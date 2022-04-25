#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject
from datetime import datetime
import logging
import os

__author__ = 'adamkoziol'


class Reporter(object):

    def metadata_reporter(self):
        """
        Creates the metadata report by pulling specific attributes from the metadata objects
        """
        logging.info('Creating summary report')
        header = '{}\n'.format(','.join(self.headers))
        # Create a string to store all the results
        data = str()
        for sample in self.metadata:
            # Add the value of the appropriate attribute to the results string
            data += GenObject.returnattr(sample, 'name')
            # SampleName
            data += GenObject.returnattr(sample.run, 'SamplePlate')
            # Genus
            data += GenObject.returnattr(sample.general, 'closestrefseqgenus')
            # SamplePurity
            data += GenObject.returnattr(sample.confindr, 'num_contaminated_snvs')
            # N50
            n50 = GenObject.returnattr(sample.quast, 'N50',
                                       number=True)
            if n50 != '-,':
                data += n50
            else:
                data += '0,'
            # NumContigs
            data += GenObject.returnattr(sample.quast, 'num_contigs',
                                         number=True)
            # TotalLength
            data += GenObject.returnattr(sample.quast, 'Total_length',
                                         number=True)
            # MeanInsertSize
            data += GenObject.returnattr(sample.quast, 'mean_insert',
                                         number=True)
            # InsertSizeSTD
            data += GenObject.returnattr(sample.quast, 'std_insert',
                                         number=True)
            # AverageCoverageDepth
            data += GenObject.returnattr(sample.qualimap, 'MeanCoveragedata',
                                         number=True)
            # CoverageDepthSTD
            data += GenObject.returnattr(sample.qualimap, 'StdCoveragedata',
                                         number=True)
            # PercentGC
            data += GenObject.returnattr(sample.quast, 'GC',
                                         number=True)
            # MASH_ReferenceGenome
            data += GenObject.returnattr(sample.mash, 'closestrefseq')
            # MASH_NumMatchingHashes
            data += GenObject.returnattr(sample.mash, 'nummatches')
            # 16S_result
            data += GenObject.returnattr(sample.sixteens_full, 'sixteens_match')
            # 16S PercentID
            data += GenObject.returnattr(sample.sixteens_full, 'percent_id')
            # CoreGenesPresent
            data += GenObject.returnattr(sample.gdcs, 'coreresults')
            # rMLST_Result
            try:
                # If the number of matches to the closest reference profile is 53, return the profile number
                if sample.rmlst.matches == 53:
                    if type(sample.rmlst.sequencetype) is list:
                        rmlst_seq_type = ';'.join(sorted(sample.rmlst.sequencetype)).rstrip(';') + ','
                    else:
                        rmlst_seq_type = GenObject.returnattr(sample.rmlst, 'sequencetype')
                        rmlst_seq_type = rmlst_seq_type if rmlst_seq_type != 'ND,' else 'new,'
                    data += rmlst_seq_type
                else:
                    # Otherwise the profile is set to new
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # MLST_Result
            try:
                if sample.mlst.matches == 7:
                    if type(sample.mlst.sequencetype) is list:
                        mlst_seq_type = ';'.join(sorted(sample.mlst.sequencetype)).rstrip(';') + ','
                    else:
                        mlst_seq_type = GenObject.returnattr(sample.mlst, 'sequencetype')
                        mlst_seq_type = mlst_seq_type if mlst_seq_type != 'ND,' else 'new,'
                    data += mlst_seq_type
                else:
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # MLST_gene_X_alleles
            try:
                # Create a set of all the genes present in the results (gene name split from allele)
                gene_set = {gene.split('_')[0] for gene in sample.mlst.combined_metadata_results}
                for gene in sorted(gene_set):
                    allele_list = list()
                    # Determine all the alleles that are present for each gene
                    for allele in sample.mlst.combined_metadata_results:
                        if gene in allele:
                            allele_list.append(allele.replace(' ', '_'))
                    # If there is more than one allele in the sample, add both to the string separated by a ';'
                    if len(allele_list) > 1:
                        data += '{},'.format(';'.join(allele_list))
                    # Otherwise add the only allele
                    else:
                        data += allele_list[0] + ','
                # If there are fewer than seven matching alleles, add a ND for each missing result
                if len(gene_set) < 7:
                    data += (7 - len(gene_set)) * 'ND,'
            except AttributeError:
                # data += '-,-,-,-,-,-,-,'
                data += 'ND,ND,ND,ND,ND,ND,ND,'
            # E_coli_Serotype
            try:
                # If no O-type was found, set the output to be O-untypeable
                if ';'.join(sample.ectyper.o_type) == '-':
                    otype = 'O-untypeable'
                else:
                    otype = sample.ectyper.o_type
                # Same as above for the H-type
                if ';'.join(sample.ectyper.h_type) == '-':
                    htype = 'H-untypeable'

                else:
                    htype = sample.ectyper.h_type
                serotype = '{otype}:{htype},'.format(otype=otype,
                                                     htype=htype)
                # Add the serotype to the data string unless neither O-type not H-type were found; add ND instead
                data += serotype if serotype != 'O-untypeable:H-untypeable,' else 'ND,'
            except AttributeError:
                data += 'ND,'
            # SISTR_serovar_antigen
            data += GenObject.returnattr(sample.sistr, 'serovar_antigen').rstrip(';')
            # SISTR_serovar_cgMLST
            data += GenObject.returnattr(sample.sistr, 'serovar_cgmlst')
            # SISTR_serogroup
            data += GenObject.returnattr(sample.sistr, 'serogroup')
            # SISTR_h1
            data += GenObject.returnattr(sample.sistr, 'h1').rstrip(';')
            # SISTR_h2
            data += GenObject.returnattr(sample.sistr, 'h2').rstrip(';')
            # SISTR_serovar
            data += GenObject.returnattr(sample.sistr, 'serovar')
            # GeneSeekr_Profile
            try:
                if sample.genesippr.report_output:
                    data += ';'.join(sample.genesippr.report_output) + ','
                else:
                    data += 'ND,'
            except AttributeError:
                data += 'ND,'
            # Vtyper_Profile
            data += GenObject.returnattr(sample.verotoxin, 'verotoxin_subtypes_set')
            # AMR_Profile and resistant/sensitive status
            if sample.resfinder_assembled.pipelineresults:
                # Profile
                for resistance, resistance_set in sorted(sample.resfinder_assembled.pipelineresults.items()):
                    data += '{res}({r_set});'.format(res=resistance.replace(',', ';'),
                                                     r_set=';'.join(sorted(list(resistance_set))))
                data += ','
                # Resistant/Sensitive
                data += 'Resistant,'
            else:
                # Profile
                data += 'ND,'
                # Resistant/Sensitive
                data += 'Sensitive,'
            # Plasmid Result'
            if sample.mobrecon.pipelineresults:
                for plasmid, details in sorted(sample.mobrecon.pipelineresults.items()):
                    data += '{plasmid}({details});'.format(plasmid=plasmid,
                                                           details=details)
                data += ','
            else:
                data += 'ND,'
            # TotalPredictedGenes
            data += GenObject.returnattr(sample.prodigal, 'predictedgenestotal',
                                         number=True)
            # PredictedGenesOver3000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover3000bp',
                                         number=True)
            # PredictedGenesOver1000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover1000bp',
                                         number=True)
            # PredictedGenesOver500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover500bp',
                                         number=True)
            # PredictedGenesUnder500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesunder500bp',
                                         number=True)
            # AssemblyDate
            data += datetime.now().strftime('%Y-%m-%d') + ','
            # PipelineVersion
            data += self.commit + ','
            # Name of the database used in the analyses
            data += os.path.split(self.reffilepath)[-1] + ','
            # Database download date
            data += self.download_date
            # Append a new line to the end of the results for this sample
            data += '\n'
        # Replace any NA values with ND
        cleandata = data.replace('NA', 'ND')
        with open(os.path.join(self.reportpath, 'combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(header)
            metadatareport.write(cleandata)

    def legacy_reporter(self):
        """
        Creates an output that is compatible with the strain database. This method will be removed once
        a new database scheme is implemented
        """
        logging.info('Creating database-friendly summary report')
        header = '{}\n'.format(','.join(self.legacy_headers))
        # Create a string to store all the results
        data = str()
        for sample in self.metadata:
            # Add the value of the appropriate attribute to the results string
            data += GenObject.returnattr(sample, 'name')
            # SampleName
            data += GenObject.returnattr(sample.run, 'SamplePlate')
            # Genus
            data += GenObject.returnattr(sample.general, 'closestrefseqgenus')
            # SequencingDate
            data += GenObject.returnattr(sample.run, 'Date')
            # Analyst
            data += GenObject.returnattr(sample.run, 'InvestigatorName')
            # Legacy ConFindr clean/contaminated call
            data += 'ND,'
            # N50
            n50 = GenObject.returnattr(sample.quast, 'N50',
                                       number=True)
            if n50 != '-,':
                data += n50
            else:
                data += '0,'
            # NumContigs
            data += GenObject.returnattr(sample.quast, 'num_contigs',
                                         number=True)
            # TotalLength
            data += GenObject.returnattr(sample.quast, 'Total_length',
                                         number=True)
            # MeanInsertSize
            data += GenObject.returnattr(sample.quast, 'mean_insert',
                                         number=True)
            # InsertSizeSTD
            data += GenObject.returnattr(sample.quast, 'std_insert',
                                         number=True)
            # AverageCoverageDepth
            data += GenObject.returnattr(sample.qualimap, 'MeanCoveragedata',
                                         number=True)
            # CoverageDepthSTD
            data += GenObject.returnattr(sample.qualimap, 'StdCoveragedata',
                                         number=True)
            # PercentGC
            data += GenObject.returnattr(sample.quast, 'GC',
                                         number=True)
            # MASH_ReferenceGenome
            data += GenObject.returnattr(sample.mash, 'closestrefseq')
            # MASH_NumMatchingHashes
            data += GenObject.returnattr(sample.mash, 'nummatches')
            # 16S_result
            data += GenObject.returnattr(sample.sixteens_full, 'sixteens_match')
            # 16S PercentID
            data += GenObject.returnattr(sample.sixteens_full, 'percent_id')
            # rMLST_Result
            try:
                # If the number of matches to the closest reference profile is 53, return the profile number
                if sample.rmlst.matches == 53:
                    if type(sample.rmlst.sequencetype) is list:
                        rmlst_seq_type = ';'.join(sorted(sample.rmlst.sequencetype)).rstrip(';') + ','
                    else:
                        rmlst_seq_type = GenObject.returnattr(sample.rmlst, 'sequencetype')
                        rmlst_seq_type = rmlst_seq_type if rmlst_seq_type != 'ND,' else 'new,'
                    data += rmlst_seq_type
                else:
                    # Otherwise the profile is set to new
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # MLST_Result
            try:
                if sample.mlst.matches == 7:
                    if type(sample.mlst.sequencetype) is list:
                        mlst_seq_type = ';'.join(sorted(sample.mlst.sequencetype)).rstrip(';') + ','
                    else:
                        mlst_seq_type = GenObject.returnattr(sample.mlst, 'sequencetype')
                        mlst_seq_type = mlst_seq_type if mlst_seq_type != 'ND,' else 'new,'
                    data += mlst_seq_type
                else:
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # MLST_gene_X_alleles
            try:
                # Create a set of all the genes present in the results (gene name split from allele)
                gene_set = {gene.split('_')[0] for gene in sample.mlst.combined_metadata_results}
                for gene in sorted(gene_set):
                    allele_list = list()
                    # Determine all the alleles that are present for each gene
                    for allele in sample.mlst.combined_metadata_results:
                        if gene in allele:
                            allele_list.append(allele.replace(' ', '_'))
                    # If there is more than one allele in the sample, add both to the string separated by a ';'
                    if len(allele_list) > 1:
                        data += '{},'.format(';'.join(allele_list))
                    # Otherwise add the only allele
                    else:
                        data += allele_list[0] + ','
                # If there are fewer than seven matching alleles, add a ND for each missing result
                if len(gene_set) < 7:
                    data += (7 - len(gene_set)) * 'ND,'
            except AttributeError:
                # data += '-,-,-,-,-,-,-,'
                data += 'ND,ND,ND,ND,ND,ND,ND,'
            # CoreGenesPresent
            data += GenObject.returnattr(sample.gdcs, 'coreresults')
            # E_coli_Serotype
            try:
                # If no O-type was found, set the output to be O-untypeable
                if ';'.join(sample.ectyper.o_type) == '-':
                    otype = 'O-untypeable'
                else:
                    otype = sample.ectyper.o_type
                # Same as above for the H-type
                if ';'.join(sample.ectyper.h_type) == '-':
                    htype = 'H-untypeable'

                else:
                    htype = sample.ectyper.h_type
                serotype = '{otype}:{htype},'.format(otype=otype,
                                                     htype=htype)
                # Add the serotype to the data string unless neither O-type not H-type were found; add ND instead
                data += serotype if serotype != 'O-untypeable:H-untypeable,' else 'ND,'
            except AttributeError:
                data += 'ND,'
            # SISTR_serovar_antigen
            data += GenObject.returnattr(sample.sistr, 'serovar_antigen').rstrip(';')
            # SISTR_serovar_cgMLST
            data += GenObject.returnattr(sample.sistr, 'serovar_cgmlst')
            # SISTR_serogroup
            data += GenObject.returnattr(sample.sistr, 'serogroup')
            # SISTR_h1
            data += GenObject.returnattr(sample.sistr, 'h1').rstrip(';')
            # SISTR_h2
            data += GenObject.returnattr(sample.sistr, 'h2').rstrip(';')
            # SISTR_serovar
            data += GenObject.returnattr(sample.sistr, 'serovar')
            # GeneSeekr_Profile
            try:
                if sample.genesippr.report_output:
                    data += ';'.join(sample.genesippr.report_output) + ','
                else:
                    data += 'ND,'
            except AttributeError:
                data += 'ND,'
            # Vtyper_Profile
            data += GenObject.returnattr(sample.verotoxin, 'verotoxin_subtypes_set')
            # AMR_Profile and resistant/sensitive status
            if sample.resfinder_assembled.pipelineresults:
                # Profile
                for resistance, resistance_set in sorted(sample.resfinder_assembled.pipelineresults.items()):
                    data += '{res}({r_set});'.format(res=resistance.replace(',', ';'),
                                                     r_set=';'.join(sorted(list(resistance_set))))
                data += ','
                # Resistant/Sensitive
                data += 'Resistant,'
            else:
                # Profile
                data += 'ND,'
                # Resistant/Sensitive
                data += 'Sensitive,'
            # Plasmid Result'
            if sample.mobrecon.pipelineresults:
                for plasmid, details in sorted(sample.mobrecon.pipelineresults.items()):
                    data += '{plasmid}({details});'.format(plasmid=plasmid,
                                                           details=details)
                data += ','
            else:
                data += 'ND,'
            # TotalPredictedGenes
            data += GenObject.returnattr(sample.prodigal, 'predictedgenestotal',
                                         number=True)
            # PredictedGenesOver3000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover3000bp',
                                         number=True)
            # PredictedGenesOver1000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover1000bp',
                                         number=True)
            # PredictedGenesOver500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover500bp',
                                         number=True)
            # PredictedGenesUnder500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesunder500bp',
                                         number=True)
            # NumClustersPF
            data += GenObject.returnattr(sample.run, 'NumberofClustersPF')
            # Percentage of reads mapping to PhiX control
            data += GenObject.returnattr(sample.run, 'phix_aligned')
            # Error rate calculated from PhiX control
            data += GenObject.returnattr(sample.run, 'error_rate')
            # LengthForwardRead
            data += GenObject.returnattr(sample.run, 'forwardlength',
                                         number=True)
            # LengthReverseRead
            data += GenObject.returnattr(sample.run, 'reverselength',
                                         number=True)
            # Real time strain
            data += GenObject.returnattr(sample.run, 'Description')
            # Flowcell
            data += GenObject.returnattr(sample.run, 'flowcell')
            # MachineName
            data += GenObject.returnattr(sample.run, 'instrument')
            # PipelineVersion
            data += self.commit + ','
            # AssemblyDate
            data += datetime.now().strftime('%Y-%m-%d') + ','
            # SamplePurity
            data += GenObject.returnattr(sample.confindr, 'num_contaminated_snvs')
            # cgMLST
            try:
                if type(sample.cgmlst.sequencetype) is list:
                    if sample.cgmlst.sequencetype:
                        cgmlst_seq_type = ';'.join(sorted(sample.cgmlst.sequencetype)).rstrip(';') + ','
                    else:
                        cgmlst_seq_type = 'ND,'
                else:
                    cgmlst_seq_type = GenObject.returnattr(sample.cgmlst, 'sequencetype')
                    # cgmlst_seq_type = cgmlst_seq_type if cgmlst_seq_type != 'ND,' else 'new,'
                data += cgmlst_seq_type
            except AttributeError:
                data += 'ND,'
            # Name of the database used in the analyses
            data += os.path.split(self.reffilepath)[-1] + ','
            # Database download date
            data += self.download_date
            # Append a new line to the end of the results for this sample
            data += '\n'
        # Replace any NA values with ND
        cleandata = data.replace('NA', 'ND')
        with open(os.path.join(self.reportpath, 'legacy_combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(header)
            metadatareport.write(cleandata)

    def clean_object(self):
        for sample in self.metadata:
            try:
                delattr(sample.coregenome, 'targetnames')
            except AttributeError:
                pass
            try:
                delattr(sample.coregenome, 'targets')
            except AttributeError:
                pass

    def run_quality_reporter(self):
        logging.info('Creating run quality summary report')
        run_name = os.path.split(self.path)[-1]
        data = 'RunName,SequencingDate,Analyst,ClusterDensity,PercentOverQ30,NumberofClustersPF,' \
               'PercentReadsPhiX,ErrorRate, LengthForwardRead,LengthReverseRead,Flowcell,MachineName\n'
        for sample in self.metadata:
            # RunName
            data += '{rn},'.format(rn=run_name)
            # SequencingDate
            data += GenObject.returnattr(sample.run, 'Date')
            # Analyst
            data += GenObject.returnattr(sample.run, 'InvestigatorName')
            # ClusterDensity
            data += GenObject.returnattr(sample.run, 'cluster_density')
            # Percentage of reads with Q-score over 30
            data += GenObject.returnattr(sample.run, 'over_q30')
            # NumClustersPF
            data += GenObject.returnattr(sample.run, 'NumberofClustersPF')
            # Percentage of reads mapping to PhiX control
            data += GenObject.returnattr(sample.run, 'phix_aligned')
            # Error rate calculated from PhiX control
            data += GenObject.returnattr(sample.run, 'error_rate')
            # LengthForwardRead
            data += GenObject.returnattr(sample.run, 'forwardlength',
                                         number=True)
            # LengthReverseRead
            data += GenObject.returnattr(sample.run, 'reverselength',
                                         number=True)
            # Flowcell
            data += GenObject.returnattr(sample.run, 'flowcell')
            # MachineName
            data += GenObject.returnattr(sample.run, 'instrument')
            break
        with open(os.path.join(self.reportpath, 'run_metrics_report.csv'), 'w') as run_report:
            run_report.write(data)

    def sample_quality_report(self):
        logging.info('Creating sample quality summary report')
        header = '{}\n'.format(','.join(self.quality_headers))
        # Create a string to store all the results
        data = str()
        for sample in self.metadata:
            # Add the value of the appropriate attribute to the results string
            data += GenObject.returnattr(sample, 'name')
            # SampleName
            data += GenObject.returnattr(sample.run, 'SamplePlate')
            # Genus
            data += GenObject.returnattr(sample.general, 'closestrefseqgenus')
            # SamplePurity
            data += GenObject.returnattr(sample.confindr, 'num_contaminated_snvs')
            # N50
            n50 = GenObject.returnattr(sample.quast, 'N50',
                                       number=True)
            if n50 != '-,':
                data += n50
            else:
                data += '0,'
            # NumContigs
            data += GenObject.returnattr(sample.quast, 'num_contigs',
                                         number=True)
            # TotalLength
            data += GenObject.returnattr(sample.quast, 'Total_length',
                                         number=True)
            # MeanInsertSize
            data += GenObject.returnattr(sample.quast, 'mean_insert',
                                         number=True)
            # InsertSizeSTD
            data += GenObject.returnattr(sample.quast, 'std_insert',
                                         number=True)
            # AverageCoverageDepth
            data += GenObject.returnattr(sample.qualimap, 'MeanCoveragedata',
                                         number=True)
            # CoverageDepthSTD
            data += GenObject.returnattr(sample.qualimap, 'StdCoveragedata',
                                         number=True)
            # PercentGC
            data += GenObject.returnattr(sample.quast, 'GC',
                                         number=True)
            # MASH_ReferenceGenome
            data += GenObject.returnattr(sample.mash, 'closestrefseq')
            # MASH_NumMatchingHashes
            data += GenObject.returnattr(sample.mash, 'nummatches')
            # rMLST_Result
            try:
                # If the number of matches to the closest reference profile is 53, return the profile number
                if sample.rmlst.matches == 53:
                    if type(sample.rmlst.sequencetype) is list:
                        rmlst_seq_type = ';'.join(sorted(sample.rmlst.sequencetype)).rstrip(';') + ','
                    else:
                        rmlst_seq_type = GenObject.returnattr(sample.rmlst, 'sequencetype')
                        rmlst_seq_type = rmlst_seq_type if rmlst_seq_type != 'ND,' else 'new,'
                    data += rmlst_seq_type
                else:
                    # Otherwise the profile is set to new
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # TotalPredictedGenes
            data += GenObject.returnattr(sample.prodigal, 'predictedgenestotal',
                                         number=True)
            # PredictedGenesOver3000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover3000bp',
                                         number=True)
            # PredictedGenesOver1000bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover1000bp',
                                         number=True)
            # PredictedGenesOver500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesover500bp',
                                         number=True)
            # PredictedGenesUnder500bp
            data += GenObject.returnattr(sample.prodigal, 'predictedgenesunder500bp',
                                         number=True)
            # PipelineVersion
            data += self.commit + ','
            # AssemblyDate
            data += datetime.now().strftime('%Y-%m-%d') + ','
            # Name of the database used in the analyses
            data += os.path.split(self.reffilepath)[-1] + ','
            # Database download date
            data += self.download_date
            # Append a new line to the end of the results for this sample
            data += '\n'
        # Replace any NA values with ND
        cleandata = data.replace('NA', 'ND')
        with open(os.path.join(self.reportpath, 'preliminary_combinedMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(header)
            metadatareport.write(cleandata)

    def lab_report(self):
        logging.info('Creating a final report')
        header = '{values}\n'.format(values=','.join(self.lab_headers))
        # Create a string to store all the results
        data = str()
        for sample in self.metadata:
            # SEQID
            data += GenObject.returnattr(sample, 'name')
            # SampleName
            data += GenObject.returnattr(sample.run, 'SamplePlate')
            # Genus
            data += GenObject.returnattr(sample.general, 'closestrefseqgenus')
            # E_coli_Serotype
            try:
                # If no O-type was found, set the output to be O-untypeable
                if ';'.join(sample.ectyper.o_type) == '-':
                    otype = 'O-untypeable'
                else:
                    otype = sample.ectyper.o_type
                # Same as above for the H-type
                if ';'.join(sample.ectyper.h_type) == '-':
                    htype = 'H-untypeable'

                else:
                    htype = sample.ectyper.h_type
                serotype = '{otype}:{htype},'.format(otype=otype,
                                                     htype=htype)
                # Add the serotype to the data string unless neither O-type not H-type were found; add ND instead
                data += serotype if serotype != 'O-untypeable:H-untypeable,' else 'ND,'
            except AttributeError:
                data += 'ND,'
            # SISTR_serovar
            data += GenObject.returnattr(sample.sistr, 'serovar')
            # GeneSeekr_Profile
            try:
                if sample.genesippr.report_output:
                    data += ';'.join(sample.genesippr.report_output) + ','
                else:
                    data += 'ND,'
            except AttributeError:
                data += 'ND,'
            # Vtyper_Profile
            data += GenObject.returnattr(sample.verotoxin, 'verotoxin_subtypes_set')
            # rMLST_Result
            try:
                # If the number of matches to the closest reference profile is 53, return the profile number
                if sample.rmlst.matches == 53:
                    if type(sample.rmlst.sequencetype) is list:
                        rmlst_seq_type = ';'.join(sorted(sample.rmlst.sequencetype)).rstrip(';') + ','
                    else:
                        rmlst_seq_type = GenObject.returnattr(sample.rmlst, 'sequencetype')
                        rmlst_seq_type = rmlst_seq_type if rmlst_seq_type != 'ND,' else 'new,'
                    data += rmlst_seq_type
                else:
                    # Otherwise the profile is set to new
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # MLST_Result
            try:
                if sample.mlst.matches == 7:
                    if type(sample.mlst.sequencetype) is list:
                        mlst_seq_type = ';'.join(sorted(sample.mlst.sequencetype)).rstrip(';') + ','
                    else:
                        mlst_seq_type = GenObject.returnattr(sample.mlst, 'sequencetype')
                        mlst_seq_type = mlst_seq_type if mlst_seq_type != 'ND,' else 'new,'
                    data += mlst_seq_type
                else:
                    data += 'new,'
            except AttributeError:
                data += 'new,'
            # N50
            n50 = GenObject.returnattr(sample.quast, 'N50',
                                       number=True)
            if n50 != '-,':
                data += n50
            else:
                data += '0,'
            # NumContigs
            data += GenObject.returnattr(sample.quast, 'num_contigs',
                                         number=True)
            # TotalLength
            data += GenObject.returnattr(sample.quast, 'Total_length',
                                         number=True)
            # AverageCoverageDepth
            data += GenObject.returnattr(sample.qualimap, 'MeanCoveragedata',
                                         number=True)
            # ConFindrContamSNVs
            data += GenObject.returnattr(sample.confindr, 'num_contaminated_snvs')
            # SequencingDate
            data += GenObject.returnattr(sample.run, 'Date')
            # Analyst
            data += GenObject.returnattr(sample.run, 'InvestigatorName')
            # Flowcell
            data += GenObject.returnattr(sample.run, 'flowcell')
            # MachineName
            data += GenObject.returnattr(sample.run, 'instrument')
            # AssemblyDate
            data += datetime.now().strftime('%Y-%m-%d') + ','
            # PipelineVersion
            data += self.commit + ','
            # Database version
            data += os.path.split(self.reffilepath)[-1] + ','
            # Append a new line to the end of the results for this sample
            data += '\n'
            # Replace any NA values with ND
        cleandata = data.replace('NA', 'ND')
        with open(os.path.join(self.reportpath, 'summaryMetadata.csv'), 'w') as metadatareport:
            metadatareport.write(header)
            metadatareport.write(cleandata)

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.commit = inputobject.commit
        self.reportpath = inputobject.reportpath
        self.starttime = inputobject.starttime
        self.path = inputobject.path
        self.reffilepath = inputobject.reffilepath
        # Define the headers to be used in the metadata report
        self.quality_headers = [
            'SeqID',
            'SampleName',
            'Genus',
            'ConfindrContamSNVs',
            'N50',
            'NumContigs',
            'TotalLength',
            'MeanInsertSize',
            'InsertSizeSTD',
            'AverageCoverageDepth',
            'CoverageDepthSTD',
            'PercentGC',
            'MASH_ReferenceGenome',
            'MASH_NumMatchingHashes',
            'rMLST_Result',
            'TotalPredictedGenes',
            'PredictedGenesOver3000bp',
            'PredictedGenesOver1000bp',
            'PredictedGenesOver500bp',
            'PredictedGenesUnder500bp',
            'AssemblyDate',
            'PipelineVersion',
            'Database',
            'DatabaseDownloadDate'
        ]
        self.headers = [
            'SeqID',
            'SampleName',
            'Genus',
            'ConfindrContamSNVs',
            'N50',
            'NumContigs',
            'TotalLength',
            'MeanInsertSize',
            'InsertSizeSTD',
            'AverageCoverageDepth',
            'CoverageDepthSTD',
            'PercentGC',
            'MASH_ReferenceGenome',
            'MASH_NumMatchingHashes',
            '16S_result',
            '16S_PercentID',
            'CoreGenesPresent',
            'rMLST_Result',
            'MLST_Result',
            'MLST_gene_1_allele',
            'MLST_gene_2_allele',
            'MLST_gene_3_allele',
            'MLST_gene_4_allele',
            'MLST_gene_5_allele',
            'MLST_gene_6_allele',
            'MLST_gene_7_allele',
            'E_coli_Serotype',
            'SISTR_serovar_antigen',
            'SISTR_serovar_cgMLST',
            'SISTR_serogroup',
            'SISTR_h1',
            'SISTR_h2',
            'SISTR_serovar',
            'GeneSeekr_Profile',
            'Vtyper_Profile',
            'AMR_Profile',
            'AMR Resistant/Sensitive',
            'PlasmidProfile',
            'TotalPredictedGenes',
            'PredictedGenesOver3000bp',
            'PredictedGenesOver1000bp',
            'PredictedGenesOver500bp',
            'PredictedGenesUnder500bp',
            'AssemblyDate',
            'PipelineVersion',
            'Database',
            'DatabaseDownloadDate'
        ]
        self.legacy_headers = [
            'SeqID',
            'SampleName',
            'Genus',
            'SequencingDate',
            'Analyst',
            'SamplePurity',
            'N50',
            'NumContigs',
            'TotalLength',
            'MeanInsertSize',
            'InsertSizeSTD',
            'AverageCoverageDepth',
            'CoverageDepthSTD',
            'PercentGC',
            'MASH_ReferenceGenome',
            'MASH_NumMatchingHashes',
            '16S_result',
            '16S_PercentID',
            'rMLST_Result',
            'MLST_Result',
            'MLST_gene_1_allele',
            'MLST_gene_2_allele',
            'MLST_gene_3_allele',
            'MLST_gene_4_allele',
            'MLST_gene_5_allele',
            'MLST_gene_6_allele',
            'MLST_gene_7_allele',
            'CoreGenesPresent',
            'E_coli_Serotype',
            'SISTR_serovar_antigen',
            'SISTR_serovar_cgMLST',
            'SISTR_serogroup',
            'SISTR_h1',
            'SISTR_h2',
            'SISTR_serovar',
            'GeneSeekr_Profile',
            'Vtyper_Profile',
            'AMR_Profile',
            'AMR Resistant/Sensitive',
            'PlasmidProfile',
            'TotalPredictedGenes',
            'PredictedGenesOver3000bp',
            'PredictedGenesOver1000bp',
            'PredictedGenesOver500bp',
            'PredictedGenesUnder500bp',
            'NumClustersPF',
            'PercentReadsPhiX',
            'ErrorRate',
            'LengthForwardRead',
            'LengthReverseRead',
            'RealTimeStrain',
            'Flowcell',
            'MachineName',
            'PipelineVersion',
            'AssemblyDate',
            'ConfindrContamSNVs',
            'cgMLST_Result',
            'Database',
            'DatabaseDownloadDate'
        ]
        self.lab_headers = [
            'SeqID',
            'SampleName',
            'Genus',
            'E_coli_Serotype',
            'SISTR_serovar',
            'GeneSeekr_Profile',
            'Vtyper_Profile',
            'rMLST_Result',
            'MLST_Result',
            'N50',
            'NumContigs',
            'TotalLength',
            'AverageCoverageDepth',
            'ConfindrContamSNVs',
            'SequencingDate',
            'Analyst',
            'Flowcell',
            'MachineName',
            'AssemblyDate',
            'PipelineVersion',
            'Database',
        ]
        try:
            with open(os.path.join(self.reffilepath, 'download_date'), 'r') as download_date:
                self.download_date = download_date.readline().rstrip()
        except FileNotFoundError:
            self.download_date = 'ND'
