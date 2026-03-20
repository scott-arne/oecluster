#include <gtest/gtest.h>
#include <atomic>
#include <vector>
#include "oecluster/ThreadPool.h"
#include "oecluster/Error.h"

using namespace OECluster;

TEST(ThreadPoolTest, ConstructDefault) {
    ThreadPool pool;
    EXPECT_GT(pool.NumThreads(), 0);
}

TEST(ThreadPoolTest, ConstructExplicitThreads) {
    ThreadPool pool(4);
    EXPECT_EQ(pool.NumThreads(), 4);
}

TEST(ThreadPoolTest, ParallelForExecutesAllItems) {
    ThreadPool pool(4);
    std::atomic<size_t> count{0};
    pool.ParallelFor(0, 1000, 64, [&count](size_t begin, size_t end) {
        count += (end - begin);
    });
    EXPECT_EQ(count.load(), 1000);
}

TEST(ThreadPoolTest, ParallelForCorrectResults) {
    ThreadPool pool(4);
    std::vector<int> results(100, 0);
    pool.ParallelFor(0, 100, 10, [&results](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            results[i] = static_cast<int>(i * 2);
        }
    });
    for (size_t i = 0; i < 100; ++i) {
        EXPECT_EQ(results[i], static_cast<int>(i * 2));
    }
}

TEST(ThreadPoolTest, CancelStopsExecution) {
    ThreadPool pool(2);
    std::atomic<size_t> count{0};
    pool.ParallelFor(0, 100000, 1, [&pool, &count](size_t, size_t) {
        count++;
        if (count > 100) pool.Cancel();
    });
    EXPECT_TRUE(pool.IsCancelled());
    EXPECT_LT(count.load(), 100000);
}

TEST(ThreadPoolTest, ExceptionPropagation) {
    ThreadPool pool(2);
    EXPECT_THROW(
        pool.ParallelFor(0, 100, 10, [](size_t begin, size_t) {
            if (begin == 50) throw OECluster::OEClusterError("test error");
        }),
        OECluster::OEClusterError
    );
}
