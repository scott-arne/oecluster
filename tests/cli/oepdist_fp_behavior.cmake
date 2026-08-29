# Asserts that oepdist "fp" options reach the distance matrix, not merely the
# JSON sidecar. oepdist_fp_metadata.cmake cannot do this: the sidecar entries
# and the FingerprintOptions assignments in tools/oepdist.cpp are independent
# statements over the same variables, so deleting an assignment leaves the
# sidecar truthful and the matrix wrong. Each pair below differs in exactly one
# option and must therefore produce a different matrix.
# Run by ctest with -DOEPDIST=<binary> -DWORK_DIR=<dir>.

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(WRITE "${WORK_DIR}/mols.smi"
     "c1ccccc1 benzene\nc1ccc(O)cc1 phenol\nCCCCCCCC octane\n")
# The chirality pair needs its own input: the molecules above are achiral, so
# --use-chirality is a no-op on them and the assertion could not fail. Morgan is
# the family used here because it is the only OEFP 0.3.0 generator that
# distinguishes enantiomers.
file(WRITE "${WORK_DIR}/chiral.smi"
     "N[C@@H](C)C(=O)O ala_r\nN[C@H](C)C(=O)O ala_s\nc1ccccc1 benzene\n")

# Runs one fp invocation over ${input} and leaves the output matrix's digest in
# ${out_var}.
function(hash_fp out_var input name)
    execute_process(
        COMMAND "${OEPDIST}" fp "${WORK_DIR}/${input}"
                -o "${WORK_DIR}/${name}.npy" ${ARGN}
        RESULT_VARIABLE status
        ERROR_VARIABLE stderr_text)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "oepdist fp ${ARGN} exited ${status}: ${stderr_text}")
    endif()
    file(MD5 "${WORK_DIR}/${name}.npy" digest)
    set(${out_var} "${digest}" PARENT_SCOPE)
endfunction()

# Fails naming the option whose value stopped reaching the matrix.
function(expect_different option first second)
    if("${first}" STREQUAL "${second}")
        message(FATAL_ERROR
                "${option}: both runs produced the same matrix (${first}). The option "
                "reaches the JSON sidecar but not the FingerprintOptions struct")
    endif()
endfunction()

hash_fp(binary_digest mols.smi binary --storage binary --metric manhattan)
hash_fp(count_digest mols.smi counted --storage count --metric manhattan)
expect_different(--storage "${binary_digest}" "${count_digest}")

hash_fp(torsion_4_digest mols.smi torsion_4
        --fp-type topological_torsions --storage count --metric manhattan
        --torsion-atom-count 4)
hash_fp(torsion_5_digest mols.smi torsion_5
        --fp-type topological_torsions --storage count --metric manhattan
        --torsion-atom-count 5)
expect_different(--torsion-atom-count "${torsion_4_digest}" "${torsion_5_digest}")

hash_fp(achiral_digest chiral.smi achiral)
hash_fp(chiral_digest chiral.smi chiral --use-chirality)
expect_different(--use-chirality "${achiral_digest}" "${chiral_digest}")
