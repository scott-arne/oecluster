/**
 * @file test_murcko_scaffold.cpp
 * @brief Tests for Bemis-Murcko scaffold extraction and clustering.
 */

#include <cstring>

#include <gtest/gtest.h>

namespace OECluster::detail {
// Declared rather than included: this commit adds no public header yet. The
// declaration is also what pulls MurckoScaffold.cpp's object out of the static
// archive, so the new OEMedChem link edge is genuinely exercised -- a stub
// translation unit nobody references would link clean either way.
const char* murcko_link_probe();
}  // namespace OECluster::detail

namespace {

TEST(MurckoBuildWiringTest, LinksAgainstOEMedChem) {
    const char* name = OECluster::detail::murcko_link_probe();
    ASSERT_NE(name, nullptr);
    // The exact spelling OEMedChem returns for a region type is not part of its
    // documented contract; that the call resolves and returns a real string is.
    EXPECT_GT(std::strlen(name), 0u);
}

}  // namespace
