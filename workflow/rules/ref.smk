rule get_ref_genome:
	output:
		temp("resources/ref_genome.fasta"),
	log:
		"logs/get_ref_genome.log",
	conda:
		"../envs/curl.yaml"

	params:
		link=config["ref_genome"]["link"],
	cache: True
	shell:
		"curl {params.link} > {output} 2> {log}"

if config["use_spikeIn"]:
	rule get_spikeIn_genome:
		output:
			temp("resources/spikeIn_genome.fasta"),
		log:
			"logs/get_spikeIn_genome.log",
		conda:
			"../envs/curl.yaml"
		params:
			link=config["spikeIn_genome"]["link"],
		cache: True
		shell:
			"curl {params.link} > {output} 2> {log}"

	rule combine_genomes:
		input:
			ref="resources/ref_genome.fasta",
			spikeIn="resources/spikeIn_genome.fasta",
		output:
			temp("resources/genome.fasta"),
		log:
			"logs/combine_genomes.log",
		cache: True
		params:
			ref_name=config["ref_genome"]["name"],
			spikeIn_name=config["spikeIn_genome"]["name"],
		shell:
			"""
			sed -e 's/>/>{params.ref_name}_/' {input.ref} > {input.ref} 2> {log}
			sed -e 's/>/>{params.spikeIn_name}_/' {input.spikeIn} > {input.spikeIn} 2>> {log}
			cat {input.ref} {input.spikeIn} > {output} 2>> {log}
			rm {input.ref} {input.spikeIn} 2>> {log}
			"""
else:
		rule rename_genome:
			input:
				"resources/ref_genome.fasta",
			output:
				"resources/genome.fasta",
			log:
				"logs/rename_genome.log",
			cache: True
			shell:
				"mv {input} {output} 2> {log}"

rule bowtie2_index:
	input:
		reference="resources/genome.fasta"
	output:
		multiext(
			"resources/genome",
			".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
		),
	log:
		"logs/bowtie2_build/build.log"
	params:
		extra=""  # optional parameters
	threads: 8
	wrapper:
		"0.77.0/bio/bowtie2/build"
		
