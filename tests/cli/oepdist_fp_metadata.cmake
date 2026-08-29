# Asserts that oepdist "fp" options reach the JSON sidecar. Exercises defaults,
# overrides, and similarity. Run by ctest with -DOEPDIST=<binary> -DWORK_DIR=<dir>.

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(WRITE "${WORK_DIR}/mols.smi"
     "c1ccccc1 benzene\nc1ccc(O)cc1 phenol\nCCCCCCCC octane\n")

# Runs one fp invocation and leaves its sidecar contents in ${out_var}.
function(run_fp out_var name)
    execute_process(
        COMMAND "${OEPDIST}" fp "${WORK_DIR}/mols.smi"
                -o "${WORK_DIR}/${name}.npy" ${ARGN}
        RESULT_VARIABLE status
        ERROR_VARIABLE stderr_text)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "oepdist fp ${ARGN} exited ${status}: ${stderr_text}")
    endif()
    file(READ "${WORK_DIR}/${name}.json" contents)
    set(${out_var} "${contents}" PARENT_SCOPE)
endfunction()

# Fails with the field name when a sidecar field is missing or wrong.
function(expect_param sidecar label field expected)
    string(JSON actual ERROR_VARIABLE json_error GET "${sidecar}" params "${field}")
    if(json_error)
        message(FATAL_ERROR "${label}: params has no \"${field}\" (${json_error})")
    endif()
    if(NOT actual STREQUAL expected)
        message(FATAL_ERROR
                "${label}: params.${field} is \"${actual}\", expected \"${expected}\"")
    endif()
endfunction()

run_fp(defaults defaults)
expect_param("${defaults}" defaults fp_type morgan)
expect_param("${defaults}" defaults storage binary)
expect_param("${defaults}" defaults numbits 2048)
expect_param("${defaults}" defaults radius 2)
expect_param("${defaults}" defaults min_distance 1)
expect_param("${defaults}" defaults max_distance 30)
expect_param("${defaults}" defaults torsion_atom_count 4)
expect_param("${defaults}" defaults use_chirality OFF)
expect_param("${defaults}" defaults metric tanimoto)
expect_param("${defaults}" defaults similarity OFF)

# Move every field except similarity off its default, so a hard-coded sidecar
# field passes the defaults block and fails here.
run_fp(overrides overrides
       --fp-type topological_torsions --storage count --numbits 1024
       --radius 3 --min-distance 2 --max-distance 12
       --torsion-atom-count 5 --use-chirality --metric manhattan)
expect_param("${overrides}" overrides fp_type topological_torsions)
expect_param("${overrides}" overrides storage count)
expect_param("${overrides}" overrides numbits 1024)
expect_param("${overrides}" overrides radius 3)
expect_param("${overrides}" overrides min_distance 2)
expect_param("${overrides}" overrides max_distance 12)
expect_param("${overrides}" overrides torsion_atom_count 5)
expect_param("${overrides}" overrides use_chirality ON)
expect_param("${overrides}" overrides metric manhattan)
expect_param("${overrides}" overrides similarity OFF)

# --sim needs its own run: it cannot ride along with the overrides above,
# because --storage count rejects tanimoto and --sim rejects manhattan. Without
# this third invocation every assertion on "similarity" would expect OFF, which
# is the default -- so a hard-coded false would pass, and that is precisely the
# regression Task 4 shipped and this test exists to prevent.
run_fp(sim sim --sim)
expect_param("${sim}" sim similarity ON)
