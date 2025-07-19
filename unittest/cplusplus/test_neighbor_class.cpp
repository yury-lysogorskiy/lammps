// unit tests for neighbor list functionality

#include "library.h"

#include "info.h"
#include "utils.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>
#include <vector>

using ::testing::ContainsRegex;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {

class NeighborListsBin : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "NeighborListsBin";
        LAMMPSTest::SetUp();
    }

    virtual void atomic_system(const std::string &atom_style, const std::string &units,
                               const std::string &newton)
    {
        BEGIN_HIDE_OUTPUT();
        command("units " + units);
        command("atom_style " + atom_style);
        command("newton " + newton);
        command("lattice sc 2.0 origin 0.125 0.125 0.125");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 2 box");
        command("create_atoms 1 box");
        command("pair_style zero 3.5");
        command("pair_coeff * *");
        command("mass * 1.0");
        command("set region box type/ratio 2 0.25 32187");
        END_HIDE_OUTPUT();
    }

    std::vector<std::string> get_neigh_info(const std::string &result)
    {
        if (verbose) utils::print(result);
        auto begin = result.find("Neighbor list info");
        auto end   = result.rfind("Per MPI rank");
        if ((begin != std::string::npos) && (end != std::string::npos) && (begin < end)) {
            return utils::split_lines(result.substr(begin, end - begin));
        } else {
            return std::vector<std::string>();
        }
    }
};

class NeighborListsNsq : public NeighborListsBin {
protected:
    void SetUp() override
    {
        testbinary = "NeighborListsNsq";
        LAMMPSTest::SetUp();
    }

    void atomic_system(const std::string &atom_style, const std::string &units,
                       const std::string &newton) override
    {
        NeighborListsBin::atomic_system(atom_style, units, newton);
        BEGIN_HIDE_OUTPUT();
        command("neighbor 2.0 nsq");
        END_HIDE_OUTPUT();
    }
};

TEST_F(NeighborListsBin, one_atomic_half_list_newton)
{
    atomic_system("atomic", "real", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style lj/cut 3.5");
    command("pair_coeff * * 0.01 2.0");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 12) {
        EXPECT_THAT(neigh_info[6], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[6], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[7], ContainsRegex(".1. pair lj/cut, perpetual"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("attributes: half, newton on"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("pair build: half/bin/atomonly/newton"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("stencil: half/bin/3d"));
        EXPECT_THAT(neigh_info[11], ContainsRegex("bin: standard"));
        int idx = lammps_find_pair_neighlist(lmp, "lj/cut", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 26);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 29);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsBin, one_atomic_half_list_nonewton)
{
    atomic_system("atomic", "real", "off");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style lj/cut 3.5");
    command("pair_coeff * * 0.01 2.0");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 12) {
        EXPECT_THAT(neigh_info[6], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[6], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[7], ContainsRegex(".1. pair lj/cut, perpetual"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("attributes: half, newton off"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("pair build: half/bin/atomonly/newtoff"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("stencil: full/bin/3d"));
        EXPECT_THAT(neigh_info[11], ContainsRegex("bin: standard"));
        int idx = lammps_find_pair_neighlist(lmp, "lj/cut", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 80);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 79);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsBin, one_atomic_full)
{
    if (!lammps_config_has_package("MANYBODY")) GTEST_SKIP() << "Missing MANYBODY package for test";
    atomic_system("atomic", "metal", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style sw");
    command("pair_coeff * * Si.sw Si Si");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 12) {
        EXPECT_THAT(neigh_info[6], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[6], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[7], ContainsRegex(".1. pair sw, perpetual"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("attributes: full, newton on"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("pair build: full/bin/atomonly"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("stencil: full/bin/3d"));
        EXPECT_THAT(neigh_info[11], ContainsRegex("bin: standard"));
        int idx = lammps_find_pair_neighlist(lmp, "sw", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 92);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 92);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsBin, one_atomic_full_ghost)
{
    if (!lammps_config_has_package("MANYBODY")) GTEST_SKIP() << "Missing MANYBODY package for test";
    atomic_system("atomic", "metal", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("neigh_modify one 5000 page 200000");
    command("pair_style airebo 3.0 1 1");
    command("pair_coeff * * CH.airebo C H");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 12) {
        EXPECT_THAT(neigh_info[6], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[6], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[7], ContainsRegex(".1. pair airebo, perpetual"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("attributes: full, newton on"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("pair build: full/bin/ghost"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("stencil: full/ghost/bin/3d"));
        EXPECT_THAT(neigh_info[11], ContainsRegex("bin: standard"));
        int idx = lammps_find_pair_neighlist(lmp, "airebo", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 948);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 948);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsBin, two_atomic_half_full)
{
    if (!lammps_config_has_package("MEAM")) GTEST_SKIP() << "Missing MEAM package for test";
    atomic_system("atomic", "metal", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style meam");
    command("pair_coeff * * library.meam Al Mg NULL Al Mg");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 17) {
        EXPECT_THAT(neigh_info[6], ContainsRegex("2 neighbor lists"));
        EXPECT_THAT(neigh_info[6], ContainsRegex("perpetual/occasional/extra = 2 0 0"));
        EXPECT_THAT(neigh_info[7], ContainsRegex(".1. pair meam, perpetual"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("attributes: full, newton on"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("pair build: full/bin/atomonly"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("stencil: full/bin/3d"));
        EXPECT_THAT(neigh_info[11], ContainsRegex("bin: standard"));
        EXPECT_THAT(neigh_info[12], ContainsRegex(".2. pair meam, perpetual"));
        EXPECT_THAT(neigh_info[13], ContainsRegex("attributes: half, newton on"));
        EXPECT_THAT(neigh_info[14], ContainsRegex("pair build: halffull/newton"));
        EXPECT_THAT(neigh_info[15], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[16], ContainsRegex("bin: none"));
        int idx = lammps_find_pair_neighlist(lmp, "meam", 1, 0, 1);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 122);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 122);
        idx = lammps_find_pair_neighlist(lmp, "meam", 1, 0, 2);
        EXPECT_EQ(idx, 1);
        inum = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 61);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 61);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsNsq, one_atomic_half_list_newton)
{
    atomic_system("atomic", "real", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style lj/cut 3.5");
    command("pair_coeff * * 0.01 2.0");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 11) {
        EXPECT_THAT(neigh_info[5], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[5], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[6], ContainsRegex(".1. pair lj/cut, perpetual"));
        EXPECT_THAT(neigh_info[7], ContainsRegex("attributes: half, newton on"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("pair build: half/nsq/newton"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("bin: none"));
        int idx = lammps_find_pair_neighlist(lmp, "lj/cut", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 46);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 53);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsNsq, one_atomic_half_list_nonewton)
{
    atomic_system("atomic", "real", "off");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style lj/cut 3.5");
    command("pair_coeff * * 0.01 2.0");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 11) {
        EXPECT_THAT(neigh_info[5], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[5], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[6], ContainsRegex(".1. pair lj/cut, perpetual"));
        EXPECT_THAT(neigh_info[7], ContainsRegex("attributes: half, newton off"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("pair build: half/nsq/newtoff"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("bin: none"));
        int idx = lammps_find_pair_neighlist(lmp, "lj/cut", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 80);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 79);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsNsq, one_atomic_full)
{
    atomic_system("atomic", "metal", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style sw");
    command("pair_coeff * * Si.sw Si Si");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 12) {
        EXPECT_THAT(neigh_info[5], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[5], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[6], ContainsRegex(".1. pair sw, perpetual"));
        EXPECT_THAT(neigh_info[7], ContainsRegex("attributes: full, newton on"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("pair build: full/nsq"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("bin: none"));
        int idx = lammps_find_pair_neighlist(lmp, "sw", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 92);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 92);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsNsq, two_atomic_half_full)
{
    atomic_system("atomic", "metal", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("pair_style meam");
    command("pair_coeff * * library.meam Al Mg NULL Al Mg");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 16) {
        EXPECT_THAT(neigh_info[5], ContainsRegex("2 neighbor lists"));
        EXPECT_THAT(neigh_info[5], ContainsRegex("perpetual/occasional/extra = 2 0 0"));
        EXPECT_THAT(neigh_info[6], ContainsRegex(".1. pair meam, perpetual"));
        EXPECT_THAT(neigh_info[7], ContainsRegex("attributes: full, newton on"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("pair build: full/nsq"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("bin: none"));
        EXPECT_THAT(neigh_info[11], ContainsRegex(".2. pair meam, perpetual"));
        EXPECT_THAT(neigh_info[12], ContainsRegex("attributes: half, newton on"));
        EXPECT_THAT(neigh_info[13], ContainsRegex("pair build: halffull/newton"));
        EXPECT_THAT(neigh_info[14], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[15], ContainsRegex("bin: none"));
        int idx = lammps_find_pair_neighlist(lmp, "meam", 1, 0, 1);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 122);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 122);
        idx = lammps_find_pair_neighlist(lmp, "meam", 1, 0, 2);
        EXPECT_EQ(idx, 1);
        inum = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 61);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 61);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}

TEST_F(NeighborListsNsq, one_atomic_full_ghost)
{
    if (!lammps_config_has_package("MANYBODY")) GTEST_SKIP() << "Missing MANYBODY package for test";
    atomic_system("atomic", "metal", "on");
    BEGIN_CAPTURE_OUTPUT();
    command("neigh_modify one 5000 page 200000");
    command("pair_style airebo 3.0 1 1");
    command("pair_coeff * * CH.airebo C H");
    command("run 0 post no");
    auto neigh_info = get_neigh_info(END_CAPTURE_OUTPUT());
    if (neigh_info.size() >= 12) {
        EXPECT_THAT(neigh_info[5], ContainsRegex("1 neighbor lists"));
        EXPECT_THAT(neigh_info[5], ContainsRegex("perpetual/occasional/extra = 1 0 0"));
        EXPECT_THAT(neigh_info[6], ContainsRegex(".1. pair airebo, perpetual"));
        EXPECT_THAT(neigh_info[7], ContainsRegex("attributes: full, newton on"));
        EXPECT_THAT(neigh_info[8], ContainsRegex("pair build: full/nsq/ghost"));
        EXPECT_THAT(neigh_info[9], ContainsRegex("stencil: none"));
        EXPECT_THAT(neigh_info[10], ContainsRegex("bin: none"));
        int idx = lammps_find_pair_neighlist(lmp, "airebo", 1, 0, 0);
        EXPECT_EQ(idx, 0);
        int nlocal = lammps_extract_setting(lmp, "nlocal");
        int inum   = lammps_neighlist_num_elements(lmp, idx);
        EXPECT_EQ(nlocal, inum);
        int numneigh   = -1;
        int iatom      = -1;
        int *neighbors = nullptr;
        lammps_neighlist_element_neighbors(lmp, idx, 0, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 0);
        EXPECT_EQ(numneigh, 948);
        lammps_neighlist_element_neighbors(lmp, idx, 1, &iatom, &numneigh, &neighbors);
        EXPECT_EQ(iatom, 1);
        EXPECT_EQ(numneigh, 948);
    } else {
        GTEST_FAIL() << "No suitable neighbor list info found";
    }
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
