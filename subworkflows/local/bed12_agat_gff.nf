include { AGAT_CONVERTBED2GFF } from '../../modules/nf-core/agat/convertbed2gff/main'
include { AGAT_CONVERTSPGXF2GXF } from '../../modules/nf-core/agat/convertspgxf2gxf/main'

workflow BED12_AGAT_GFF {
    take:
    bed                     // channel: [ meta, .bed ]

    main:
    ch_versions             = Channel.empty()

    // MODULE: AGAT_CONVERTBED2GFF
    AGAT_CONVERTBED2GFF ( bed )

    ch_versions             = ch_versions.mix(AGAT_CONVERTBED2GFF.out.versions.first())

    // Clean-up Gff attributes
    ch_merged_gff           = AGAT_CONVERTBED2GFF.out.gff
                            | map { meta, gff ->
                                def gff_lines = gff.readLines()

                                def feat_parent_id = null
                                def modified_lines = gff_lines
                                    .collect { line ->
                                        if (line.startsWith('#')) {
                                            return line
                                        }

                                        def fields = line.tokenize('\t')

                                        if (fields[2] != 'mRNA') {
                                            return fields[0..7].join('\t') + '\t' + 'Parent=' + feat_parent_id
                                        }

                                        def attributes = fields[8].tokenize(';')

                                        def name_match = attributes =~ /Name=([^;,]*)/
                                        def tama_tx_id = java.net.URLDecoder.decode(name_match[0][1], "UTF-8").tokenize(';')[1].strip()
                                        def gene_id = tama_tx_id.tokenize('.')[0]
                                        feat_parent_id = tama_tx_id

                                        return fields[0..7].join('\t') + '\t' + 'ID=' + tama_tx_id + ';Parent=' + gene_id
                                    }

                                [ "${meta.id}.minimal.gff" ] + modified_lines.join('\n')
                            }
                            | collectFile(newLine: true)
                            | map { file ->
                                [ [ id: file.baseName.replace('.minimal', '') ], file ]
                            }

    // MODULE: AGAT_CONVERTSPGXF2GXF
    AGAT_CONVERTSPGXF2GXF ( ch_merged_gff )
    ch_versions             = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    emit:

    versions                = ch_versions                           // channel: [ versions.yml ]
    gff                     = AGAT_CONVERTSPGXF2GXF.out.output_gff  // channel: [ meta, .gff ]
}
