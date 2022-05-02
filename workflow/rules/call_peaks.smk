rule macs2_call_peaks_narrow:
	input:
		unpack(get_macs2_input_narrow)
	output:
		multiext("results/peaks/individual/narrow/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/callpeak_narrow_{sample}.log"
	params: 
		config["params"]["macs2_call_peaks_narrow"],
		macs2_read_format
	wrapper:
		"v1.1.0/bio/macs2/callpeak"

rule macs2_call_peaks_broad:
	input:
		unpack(get_macs2_input_broad)
	output:
		multiext("results/peaks/individual/broad/{sample}",
                 "_peaks.xls",   ### required
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
	log:
		"logs/macs2/callpeak_broad_{sample}.log"
	params: 
		config["params"]["macs2_call_peaks_broad"],
		macs2_read_format
	wrapper:
		"v1.1.0/bio/macs2/callpeak"
		
rule macs2_call_peaks_narrow_merged:
	input:
		unpack(get_macs2_input_narrow_merged)
	output:
		multiext("results/peaks/merged/narrow/{sample_group}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/callpeak_narrow_{sample_group}.log"
	params: 
		config["params"]["macs2_call_peaks_narrow"],
		macs2_read_format_merged
	wrapper:
		"v1.1.0/bio/macs2/callpeak"

rule macs2_call_peaks_broad_merged:
	input:
		unpack(get_macs2_input_broad_merged)
	output:
		multiext("results/peaks/merged/broad/{sample_group}",
                 "_peaks.xls",   ### required
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
	log:
		"logs/macs2/callpeak_broad_{sample_group}.log"
	params: 
		config["params"]["macs2_call_peaks_broad"],
		macs2_read_format_merged
	wrapper:
		"v1.1.0/bio/macs2/callpeak"

rule filter_peaks_by_replicates:
	input:
		ind_peaks=get_replicate_peaks,
		merged_peaks="results/peaks/merged/{peak_type}/{sample_group}_peaks.{peak_type}Peak",
	output:
		"results/peaks/filtered/{sample_group}.{peak_type}Peak"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	params:
		min_overlap=config["params"]["filter_peaks_by_replicates"]["min_overlap"],
		replicate_threshold=config["params"]["filter_peaks_by_replicates"]["replicate_threshold"]
	script:
		"../scripts/filter_peaks_by_replicates.R"
