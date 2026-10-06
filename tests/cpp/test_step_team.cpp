/**
 * @file test_step_team.cpp
 * @brief The persistent team behind the spanning-tree kernel.
 */

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <vector>

#include "../../src/clustering/StepTeam.h"

using OECluster::detail::StepTeam;
using OECluster::detail::resolve_participants;
using OECluster::detail::work_unit;

TEST(StepTeamTest, ResolvesParticipantsWithAOneWorkerFloor) {
    EXPECT_EQ(resolve_participants(0, 0, 8), 0u);
    // hardware_concurrency() may report 0; the run still needs one worker.
    EXPECT_EQ(resolve_participants(0, 10, 0), 1u);
    EXPECT_EQ(resolve_participants(0, 10, 8), 8u);
    EXPECT_EQ(resolve_participants(0, 3, 8), 3u);
    EXPECT_EQ(resolve_participants(5, 10, 0), 5u);
    EXPECT_EQ(resolve_participants(50, 10, 8), 10u);
}

TEST(StepTeamTest, WorkUnitsSpreadTheRangeAndRespectTheCeiling) {
    EXPECT_EQ(work_unit(0, 4), 1u);
    EXPECT_EQ(work_unit(1, 4), 1u);
    EXPECT_EQ(work_unit(160, 4), 10u);
    EXPECT_EQ(work_unit(161, 4), 11u);
    EXPECT_EQ(work_unit(160, 4, 3), 3u);
    EXPECT_EQ(work_unit(160, 0), 40u);
    const size_t big = std::numeric_limits<size_t>::max();
    EXPECT_EQ(work_unit(big, big), 1u);
    EXPECT_EQ(work_unit(big, 1, 64), 64u);
}

TEST(StepTeamTest, RunsEveryParticipantOncePerStep) {
    for (size_t participants : {1, 2, 4, 8}) {
        StepTeam team(participants);
        EXPECT_EQ(team.Participants(), participants);
        for (int step = 0; step < 200; ++step) {
            std::vector<std::atomic<int>> calls(participants);
            team.Run([&](size_t participant) { calls[participant].fetch_add(1); });
            for (size_t p = 0; p < participants; ++p) {
                ASSERT_EQ(calls[p].load(), 1) << participants << " participants, step " << step;
            }
        }
    }
}

TEST(StepTeamTest, TheCallerIsParticipantZero) {
    StepTeam team(4);
    const std::thread::id caller = std::this_thread::get_id();
    std::thread::id seen;
    team.Run([&](size_t participant) {
        if (participant == 0) {
            seen = std::this_thread::get_id();
        }
    });
    EXPECT_EQ(seen, caller);
}

TEST(StepTeamTest, AWorkerExceptionPropagatesAndTheTeamKeepsWorking) {
    StepTeam team(4);
    EXPECT_THROW(team.Run([](size_t participant) {
                     if (participant == 3) {
                         throw std::runtime_error("worker failed");
                     }
                 }),
                 std::runtime_error);
    EXPECT_FALSE(team.Failed());
    std::atomic<int> calls{0};
    team.Run([&](size_t) { calls.fetch_add(1); });
    EXPECT_EQ(calls.load(), 4);
}

TEST(StepTeamTest, FailedStopsTheOtherParticipantsEarly) {
    StepTeam team(4);
    std::atomic<bool> saw_failure{false};
    EXPECT_THROW(team.Run([&](size_t participant) {
                     if (participant == 0) {
                         throw std::runtime_error("caller failed");
                     }
                     const auto start = std::chrono::steady_clock::now();
                     while (!team.Failed()) {
                         if (std::chrono::steady_clock::now() - start >
                             std::chrono::seconds(5)) {
                             return;
                         }
                         std::this_thread::yield();
                     }
                     saw_failure.store(true);
                 }),
                 std::runtime_error);
    EXPECT_TRUE(saw_failure.load());
}

// Idle threads sleep after spinning and yielding; a step published after that
// must still wake them.
TEST(StepTeamTest, SleepingWorkersWakeForTheNextStep) {
    StepTeam team(3);
    std::atomic<int> calls{0};
    team.Run([&](size_t) { calls.fetch_add(1); });
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    team.Run([&](size_t) { calls.fetch_add(1); });
    EXPECT_EQ(calls.load(), 6);
}

TEST(StepTeamTest, ASpawnFailureStopsTheWaitingThreads) {
    std::atomic<size_t> spawned{0};
    std::atomic<size_t> exited{0};
    auto spawn = [&](std::function<void()> body) {
        if (spawned.fetch_add(1) == 2) {
            // The first two are already waiting for a step that never comes.
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
            throw std::system_error(
                std::make_error_code(std::errc::resource_unavailable_try_again));
        }
        return std::thread([body = std::move(body), &exited] {
            body();
            exited.fetch_add(1);
        });
    };
    try {
        StepTeam team(5, spawn);
        FAIL() << "Expected std::system_error";
    } catch (const std::system_error& error) {
        EXPECT_EQ(error.code(), std::errc::resource_unavailable_try_again);
    }
    EXPECT_EQ(spawned.load(), 3u);
    EXPECT_EQ(exited.load(), 2u);
}
