# Asserts that "oepdist fp" refuses an option the rest of the configuration
# would ignore, matching the rules the Python surface applies in
# reject_inapplicable_fingerprint_kwargs (python/oecluster/_comparisons.py).
# Before those rules reached the CLI, "--max-distance 4 --fp-type morgan"
# produced a matrix byte-identical to a run without the flag -- and on 4.2.3
# that same flag was the Morgan radius, so a 4.x script got a different
# fingerprint than it asked for, silently.
# Run by ctest with -DOEPDIST=<binary> -DWORK_DIR=<dir>.

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(WRITE "${WORK_DIR}/mols.smi"
     "c1ccccc1 benzene\nc1ccc(O)cc1 phenol\nCCCCCCCC octane\n")

# Fails unless the invocation exits nonzero with ${expected} on stderr. Both
# halves matter: an exit status alone would be satisfied by any unrelated
# failure, and a message alone would be satisfied by a warning.
function(expect_refused label expected)
    execute_process(
        COMMAND "${OEPDIST}" fp "${WORK_DIR}/mols.smi"
                -o "${WORK_DIR}/${label}.npy" ${ARGN}
        RESULT_VARIABLE status
        OUTPUT_VARIABLE stdout_text
        ERROR_VARIABLE stderr_text)
    if(status EQUAL 0)
        message(FATAL_ERROR "${label}: oepdist fp ${ARGN} succeeded; expected a refusal")
    endif()
    string(FIND "${stderr_text}" "${expected}" found)
    if(found EQUAL -1)
        message(FATAL_ERROR
                "${label}: stderr does not contain \"${expected}\".\n${stderr_text}")
    endif()
endfunction()

# Fails unless the invocation succeeds. Guards the other direction: a rule that
# refused its owning family too would still pass every assertion above.
function(expect_accepted label)
    execute_process(
        COMMAND "${OEPDIST}" fp "${WORK_DIR}/mols.smi"
                -o "${WORK_DIR}/${label}.npy" ${ARGN}
        RESULT_VARIABLE status
        ERROR_VARIABLE stderr_text)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "${label}: oepdist fp ${ARGN} exited ${status}: ${stderr_text}")
    endif()
endfunction()

# --- The four per-family options, each refused by a family that ignores it ---

expect_refused(radius_on_atom_pair
               "--radius does not apply to --fp-type atom_pair"
               --fp-type atom_pair --radius 3)
expect_accepted(radius_on_morgan --fp-type morgan --radius 3)

expect_refused(min_distance_on_morgan
               "--min-distance does not apply to --fp-type morgan"
               --fp-type morgan --min-distance 2)
expect_accepted(min_distance_on_atom_pair --fp-type atom_pair --min-distance 2)

expect_refused(max_distance_on_morgan
               "--max-distance does not apply to --fp-type morgan"
               --fp-type morgan --max-distance 4)
expect_accepted(max_distance_on_atom_pair --fp-type atom_pair --max-distance 4)

expect_refused(torsion_atom_count_on_morgan
               "--torsion-atom-count does not apply to --fp-type morgan"
               --fp-type morgan --torsion-atom-count 5)
expect_accepted(torsion_atom_count_on_torsions
                --fp-type topological_torsions --torsion-atom-count 5)

# The remedy has to name a family that actually reads the option. The message
# for --max-distance under morgan offers --radius, which morgan does read, and
# --fp-type atom_pair, which reads --max-distance.
expect_refused(remedy_names_a_single_family
               "Use --radius instead, or select --fp-type atom_pair."
               --fp-type morgan --max-distance 4)

# --- The default value is not an absent option ---

# 30 is the --max-distance default, so a rule that compared values instead of
# counting occurrences would let this through. The caller asked for an
# atom-pair window on a family that has none either way.
expect_refused(explicit_default_still_refused
               "--max-distance does not apply to --fp-type morgan"
               --fp-type morgan --max-distance 30)

# --- Aliases fold onto the canonical family before the rule fires ---

# atompair and topological_atom_pair are the same generator as atom_pair; a
# rule keyed on the raw spelling would refuse --max-distance here.
expect_accepted(max_distance_on_atompair_alias --fp-type atompair --max-distance 4)
expect_accepted(max_distance_on_topological_alias
                --fp-type topological_atom_pair --max-distance 4)
expect_refused(radius_on_atompair_alias
               "--radius does not apply to --fp-type atompair"
               --fp-type atompair --radius 3)

# --- numbits against a storage with no width to fold to ---

expect_refused(numbits_on_sparse
               "--numbits does not apply to --storage sparse"
               --storage sparse --numbits 1024)
expect_refused(numbits_on_sparse_count
               "--numbits does not apply to --storage sparse_count"
               --storage sparse_count --numbits 1024 --metric manhattan)
expect_accepted(numbits_on_binary --storage binary --numbits 1024)
expect_accepted(numbits_on_count --storage count --numbits 1024 --metric manhattan)
# Sparse storage is fine on its own; it is only the width that has no meaning.
expect_accepted(sparse_without_numbits --storage sparse)

# --- An unrecognized value gets the library's message, not an advisory one ---

# Both rules stand aside on a family or storage they cannot name, so the
# constructor speaks first about the thing that is actually wrong. Without
# that, "--fp-type morgna --radius 3" would be refused for the radius.
expect_refused(unknown_family_outranks_the_radius_rule
               "Unknown OEFP fingerprint type"
               --fp-type morgna --radius 3)
expect_refused(unknown_storage_outranks_the_numbits_rule
               "Unknown fingerprint storage"
               --storage sparse_binary --numbits 1024)

# --- The refusal comes before any work ---

# The rules read only what was typed, so nothing should have been written.
if(EXISTS "${WORK_DIR}/max_distance_on_morgan.npy")
    message(FATAL_ERROR
            "max_distance_on_morgan: a matrix was written for a refused invocation")
endif()
