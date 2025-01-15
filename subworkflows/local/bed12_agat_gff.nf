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

    // MODULE: AGAT_CONVERTSPGXF2GXF
    AGAT_CONVERTSPGXF2GXF ( AGAT_CONVERTBED2GFF.out.gff )
    ch_versions             = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    emit:

    versions                = ch_versions                           // channel: [ versions.yml ]
    gff                     = AGAT_CONVERTSPGXF2GXF.out.output_gff  // channel: [ meta, .gff ]
}
