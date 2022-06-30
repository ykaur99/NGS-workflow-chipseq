import pandas as pd


# def get_genome_fn():
# 	if config["use_spikeIn"]:
# 		return "resources/genome.fasta"
# 	else:
# 		return "resources/ref_genome.fasta"
# 		
# def get_genome_bn():
# 	if config["use_spikeIn"]:
# 		return "resources/genome"
# 	else:
# 		return "resources/ref_genome"

# read in table with sample metadata
# samples = (
#     pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
#     .set_index("sample_name", drop=False)
#     .sort_index()
# )



units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

# function to check config files for inclusion of optional workflow steps
def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))

def get_fq_merge(wildcards):
	unit = units.loc[wildcards.sample]
	if all(pd.isna(unit["fq1"])):
		accession = unit["sra"]
		if all(unit["read_format"] == "SE"):
			return expand(
				"data/sra/se/{accession}.fastq.gz", accession=accession)
		else:
			return expand(
				"data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read = wildcards.read[-1])
	if all(unit["read_format"] == "SE"):
		return units.loc[wildcards.sample, "fq1"].tolist()
	fq = "fq{}".format(wildcards.read[-1])
	return units.loc[wildcards.sample, fq].tolist()


def get_bowtie2_input(wildcards):
	if not is_activated("mergeReads"):
		unit = units.loc[wildcards.sample]
		if all(pd.isna(unit["fq1"])):
			# SRA sample (always paired-end for now)
			accession = unit["sra"]
			if all(unit["read_format"] == "SE"):
				return expand("data/sra/se/{accession}.fastq.gz", accession=accession)
			else:
				return expand("data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read=[1,2])
		fastqs = units.loc[wildcards.sample, ["fq1", "fq2"]].squeeze().dropna()
		if len(fastqs) == 2:
			return [fastqs.fq1, fastqs.fq2]
		return [fastqs.fq1]
	unit = units.loc[wildcards.sample]
	if all(unit["read_format"] == "SE"):
		return ["data/merged/{sample}_single.fastq.gz"]
	return ["data/merged/{sample}_1.fastq.gz", "data/merged/{sample}_2.fastq.gz"]

def get_bam_merge(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	group = pd.unique(unit["sample_name"])
	return expand(
		"results/aligned_reads/filtered/{group}.bam", group=group)

def macs2_read_format(wildcards):
		unit = units.loc[wildcards.sample]
		if all(unit["call_peaks"]):
			if all(unit["read_format"] == "SE"):
				return "-f BAM"
			elif all(unit["read_format"] == "PE"):
				return "-f BAMPE"

def macs2_read_format_merged(wildcards):
		unit =  units[units["sample_group"] == wildcards.sample_group]
		if all(unit["call_peaks"]):
			if all(unit["read_format"] == "SE"):
				return "-f BAM"
			elif all(unit["read_format"] == "PE"):
				return "-f BAMPE"


def get_macs2_input_narrow(wildcards):
	unit = units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		if all(unit["peak_type"] == "narrow"):
			if all(pd.isna(unit["input"])):
				return {"treatment": "results/aligned_reads/filtered/{sample}.bam"}
			else:
				return {"treatment": "results/aligned_reads/filtered/{sample}.bam", "control": "results/aligned_reads/filtered/{input}.bam".format(input=unit.iloc[0].input)}

def get_macs2_input_broad(wildcards):
	unit = units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		if all(unit["peak_type"] == "broad"):
			if all(pd.isna(unit["input"])):
				return {"treatment": "results/aligned_reads/filtered/{sample}.bam"}
			else:
				return {"treatment": "results/aligned_reads/filtered/{sample}.bam", "control": "results/aligned_reads/filtered/{input}.bam".format(input=unit.iloc[0].input)}

def get_macs2_input_narrow_merged(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	if all(unit["call_peaks"]):
		if all(unit["peak_type"] == "narrow"):
			if all(pd.isna(unit["input"])):
				sample_names = pd.unique(unit["sample_name"])
				return {"treatment": expand("results/aligned_reads/filtered/{sample}.bam",sample=sample_names)}
			else:
				sample_names = pd.unique(unit["sample_name"])
				input_names = pd.unique(units.loc[unit["input"]]["sample_name"])
				return {"treatment": expand("results/aligned_reads/filtered/{sample}.bam",sample=sample_names), "control": expand("results/aligned_reads/filtered/{input}.bam",input=input_names)}

def get_macs2_input_broad_merged(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	if all(unit["call_peaks"]):
		if all(unit["peak_type"] == "broad"):
			if all(pd.isna(unit["input"])):
				sample_names = pd.unique(unit["sample_name"])
				return {"treatment": expand("results/aligned_reads/filtered/{sample}.bam",sample=sample_names)}
			else:
				sample_names = pd.unique(unit["sample_name"])
				input_names = pd.unique(units.loc[unit["input"]]["sample_name"])
				return {"treatment": expand("results/aligned_reads/filtered/{sample}.bam",sample=sample_names), "control": expand("results/aligned_reads/filtered/{input}.bam",input=input_names)}

def get_replicate_peaks(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	group = pd.unique(unit["sample_name"])
	return expand(
		"results/peaks/individual/{peak_type}/{group}_peaks.{peak_type}Peak", group=group, peak_type=wildcards.peak_type)
	


def get_scaling_input(wildcards):
	stat_files = expand(
				["results/aligned_reads/stats/{sample}_unireads.idxstats"],
				sample = units["sample_name"]
			)
	return stat_files

def get_ind_spikeIn_input(wildcards):
	unit=units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		return "results/bigwigs/zscore_normalized/individual/{sample}.bw".format(sample = wildcards.sample)

def get_merged_spikeIn_input(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample]
	if all(unit["call_peaks"]):
		return "results/bigwigs/zscore_normalized/merged/{sample}.bw".format(sample = wildcards.sample)


def get_final_output():
	final_output = []

	
		# z-score normalized bigwigs for individual replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/individual/{sample}.bw"
					],
					sample = units["sample_name"]
				)
			)

	
	# z-score normalized bigwigs for merged replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/merged/{sample}.bw"
					],
					sample = units["sample_group"]
				)
			)

	if config["use_spikeIn"]:
		# spikeIn-normalized bigwigs for individual replicates
		final_output.extend(expand(
						[
							"results/bigwigs/spikeIn_normalized/individual/{sample}.bw"
						],
						sample = units.loc[units["call_peaks"],"sample_name"]
					)
				)

		# spikeIn-normalized bigwigs for merged replicates
		final_output.extend(expand(
						[
							"results/bigwigs/spikeIn_normalized/merged/{sample}.bw"
						],
						sample = units.loc[units["call_peaks"],"sample_group"]
					)
				)



	# add narrow peak output
	if any( (units["call_peaks"]) & (units["peak_type"] == "narrow")):
		out_samples =  units[(units["call_peaks"]) & (units["peak_type"] == "narrow")]
		final_output.extend(expand(
				[
					"results/peaks/individual/narrow/{sample}{ext}"
				],
				sample = out_samples["sample_name"],
				ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
			)
		)
	
	if any( (units["call_peaks"]) & (units["peak_type"] == "narrow")):
		out_samples =  units[(units["call_peaks"]) & (units["peak_type"] == "narrow")]
		replicate_summary = out_samples.groupby('sample_group').count()
		samples_w_replicates = replicate_summary
		if samples_w_replicates.shape[0] > 0:
			out_sample_groups = samples_w_replicates.index.values.tolist()
			final_output.extend(expand(
					[
						"results/peaks/merged/narrow/{sample}{ext}"
					],
					sample = out_sample_groups,
					ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
				)
			)
	
	# add broad peak output
	if any( (units["call_peaks"]) & (units["peak_type"] == "broad")):
		out_samples =  units[(units["call_peaks"]) & (units["peak_type"] == "broad")]
		final_output.extend(expand(
				[
					"results/peaks/individual/broad/{sample}{ext}"
				],
				sample = out_samples["sample_name"],
				ext = ["_peaks.xls", "_peaks.broadPeak","_peaks.gappedPeak"]
			)
		)

	if any( (units["call_peaks"]) & (units["peak_type"] == "broad")):
		out_samples =  units[(units["call_peaks"]) & (units["peak_type"] == "broad")]
		replicate_summary = out_samples.groupby('sample_group').count()
		samples_w_replicates = replicate_summary
		if samples_w_replicates.shape[0] > 0:
			out_sample_groups = samples_w_replicates.index.values.tolist()
			final_output.extend(expand(
					[
						"results/peaks/merged/broad/{sample}{ext}"
					],
					sample = out_sample_groups,
					ext = ["_peaks.xls", "_peaks.broadPeak","_peaks.gappedPeak"]
				)
			)
	
	# add filtered peak output
	if any( (units["call_peaks"]) & (units["peak_type"] == "narrow")):
		out_samples =  units[(units["call_peaks"]) & (units["peak_type"] == "narrow")]
		replicate_summary = out_samples.groupby('sample_group').count()
		samples_w_replicates = replicate_summary
		if samples_w_replicates.shape[0] > 0:
			out_sample_groups = samples_w_replicates.index.values.tolist()
			final_output.extend(expand(
					[
						"results/peaks/filtered/{sample}.narrowPeak"
					],
					sample = out_sample_groups
				)
			)
	
	if any( (units["call_peaks"]) & (units["peak_type"] == "broad")):
		out_samples =  units[(units["call_peaks"]) & (units["peak_type"] == "broad")]
		replicate_summary = out_samples.groupby('sample_group').count()
		samples_w_replicates = replicate_summary
		if samples_w_replicates.shape[0] > 0:
			out_sample_groups = samples_w_replicates.index.values.tolist()
			final_output.extend(expand(
					[
						"results/peaks/filtered/{sample}.broadPeak"
					],
					sample = out_sample_groups
				)
			)

	
	return final_output