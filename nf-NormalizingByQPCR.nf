/*  nohup nextflow run -c nf-NormalizingByQPCR.config --name 'Dummy' nf-NormalizingByQPCR.nf -with-trace -resume -bg --input_design "DesignQpcr.csv"

*/
Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName,  
                    file("$params.input_dir/$row.BamFile", checkIfExists: true),
                    file("$params.input_dir/${row.BamFile}.bai", checkIfExists: true),
                    file("$params.input_dir/$row.BwFile", checkIfExists: true), 
                    file("$params.input_dir/$row.QpcrFile", checkIfExists: true), 
                    "$row.SampleName", "$row.CtrlName"] }
      .into { design_qpcr_csv;test1_ch }

//test1_ch.view()

/* TO DO 
VERIFY qPCR.tsv files and make sure the names match with the ones in the bed file.
publishDir : 
- bw normalise
- amplicon values
- table of normalization
- normalization factor + description
*/

process Amplicon_TagDensity {
    tag "$LibName"
    publishDir "${params.outdir}/TagDensity", mode: 'copy', //params.publish_dir_mode,
    saveAs: { filename ->
            if (filename.endsWith('.amplicon.bed')) "./$filename"
            else null
    }

    input:
    tuple LibName, file(BamFile),file(BaiFile), file(BwFile), file(QpcrFile), SampleName, CtrlName from design_qpcr_csv
    path pcr_regions from params.amplicon_pcr


    output:
    tuple LibName, file(BamFile),file(BaiFile), file(BwFile), file(QpcrFile), file("${LibName}.amplicon.bed"), SampleName, CtrlName into (density_pcr_ngs, control_ch, sample_ch)
    file(temp_file)

    script:
    """
    get_tag_density -f ${BwFile} ${pcr_regions} > temp_file
    awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$8}' temp_file > ${LibName}.amplicon.bed
    """
}
control_ch
    .filter { it[7] == "" }
    .map { it -> [ it[6],it[0],it[4],it[5] ]} //Only using LibName, QpcrFile and NgsFile
    .into { control_qpcr_ch; test2 }
//test2.view()

sample_ch
    .filter { it[7] != "" }
    .map { it -> [ it[7],it[0],it[1],it[2],it[3],it[4],it[5] ]} // keeping all of it
    .into { sample_qpcr_ch; test3 }
//test3.view()

/*Crossing control channel with sample channel
    using the CtrlName-SampleName as key to cross the two channels
    re-mapping the channel to have only one tuple instead of two.
    This re-mapping is mandatory(?) to be able to declare "file()" in the following process. 
    This declaration is required for nextflow to establish symbolic links.
    Since we are using a singularity container, links must be established in the working directory so that singularity has the reading rights.
*/
control_qpcr_ch
 .cross(sample_qpcr_ch) // USE PARENTHESIS IN CROSS FUNCTION
 .map { it-> [ it[0][0], it[0][1],it[0][2],it[0][3],  it[1][1], it[1][2], it[1][3], it[1][4], it[1][5] , it[1][6] ]} // rewriting the channel to have just one array instead of 2
 .into{ toprocess_qpcr_ch; test4}

//test4.view()


/*R Script that compares Sample/Ref_NGS & Sample/Ref_QPCR
    -import files.
    -remove non alphanumeric symboles from names
    -get the list of regions for which data is available in all files
    -compute Sample/Ctrl ratios for both qPCR and NGS
    - Outputs csv file 
    - Outputs the Normalization factor for Sample_NGS as STDOUT
Add the Normalization factor to the channel*/

process R_qpcrNorm {
    tag "$LibName"
    publishDir "${params.outdir}/NormalizationTable", mode: 'copy', //params.publish_dir_mode,
    saveAs: { filename ->
            if (filename.endsWith('.Normalization_table.tsv')) "./$filename"
            else null
    }

    container ='' //set to nothing or to a container containing R.
    input:
    tuple CtrlName, CtrlLibName, file(CtrlQPCR_file), file(CtrlNGS_File), LibName, file(BamFile), file(BaiFile), file(BwFile), file(QpcrFile),file(NgsFile) from toprocess_qpcr_ch
    output:
    tuple LibName, file(BamFile),file(BaiFile),  stdout into tonormalize_bw_ch
    file("${LibName}.Normalization_table.tsv")
    file("r_file_2_run.R")
    script:
    """
    echo "#!/usr/bin/env Rscript
    data=list(  
        CtrlQpcr=read.table('${CtrlQPCR_file}', header=FALSE, stringsAsFactor=FALSE),
        CtrlNGS=read.table('${CtrlNGS_File}', header=FALSE, stringsAsFactor=FALSE)[,c('V4', 'V7')],
        SampleQpcr=read.table('${QpcrFile}', header=FALSE, stringsAsFactor=FALSE),
        SampleNGS=read.table('${NgsFile}', header=FALSE, stringsAsFactor=FALSE)[,c('V4', 'V7')])
    original_names=data[[3]][,1]
    for(i in 1:4){colnames(data[[i]])=c('region', 'value');data[[i]][,1]=gsub(data[[i]][,1], pattern='\\\\\\\\W', perl=TRUE, replacement='')}
    original_names=cbind(original_names, names(data[['SampleQpcr']][,1]))
    Qpcr_data=merge(x=data[['SampleQpcr']], y=data[['CtrlQpcr']], by='region', all=FALSE, suffixes=c('_Qpcr_Sample', '_Qpcr_Ctrl'))
    Qpcr_data[['ratio_Qpcr']]=Qpcr_data[,2]/Qpcr_data[,3]
    Ngs_data=merge(x=data[['SampleNGS']], y=data[['CtrlNGS']], by='region', all=FALSE, suffixes=c('_Ngs_Sample', '_Ngs_Ctrl'))
    Ngs_data[['ratio_NGS']]=Ngs_data[,2]/Ngs_data[,3]
    Ngs_data=Ngs_data[ Ngs_data[,'region'] %in% Qpcr_data[,'region'], ]
    Qpcr_data=merge(Qpcr_data, Ngs_data, by='region')
    write.table(Qpcr_data, file='${LibName}.Normalization_table.tsv', sep='\\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
    factor=mean(Qpcr_data[['ratio_NGS']]/Qpcr_data[['ratio_Qpcr']])
    cat(factor, file=stdout())" > r_file_2_run.R
    Rscript --vanilla r_file_2_run.R
    """
}

/*Process that creates bigwig file from BW & Normalization factor.*/
process Normalized_bamCoverage {
    tag "$LibName"
    label "multiCpu"
    publishDir "${params.outdir}/GenomeCoverage", mode: 'copy', //params.publish_dir_mode,
    saveAs: { filename ->
            if (filename.endsWith('.bw')) "./$filename"
            else null
    }
    input: 
    tuple LibName, file(BamFile),file(BaiFile), NormFactor from tonormalize_bw_ch
    output:
    file("${LibName}.${params.mapper_id}.${params.genome_prefix}.bin${params.bin_size}.RPM.QPCR.bamCoverage.bw")
   """
   bamCoverage \
   -b ${BamFile} \
   -o ${LibName}.${params.mapper_id}.${params.genome_prefix}.bin${params.bin_size}.RPM.QPCR.bamCoverage.bw -of bigwig \
   --extendReads --centerReads --samFlagInclude 3 \
   --scaleFactor ${NormFactor} --normalizeUsing CPM  --exactScaling \
   --binSize ${params.bin_size} -p ${task.cpus}
   """
}
