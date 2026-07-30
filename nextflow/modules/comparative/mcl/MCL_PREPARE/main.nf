process MCL_PREPARE {
    label    'comparative_mcl'
    tag      "mcl_prepare"

    input:
    path blastp_tsv

    output:
    path "mcl_input.abc", emit: abc

    script:
    """
    awk -F'\\t' '\$1 != \$2 {print \$1, \$2, \$12}' ${blastp_tsv} > mcl_input.abc
    """

    stub:
    """
    printf 'A B 200\n' > mcl_input.abc
    """
}
