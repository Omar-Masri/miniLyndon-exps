reference = config["reference"]
reads, exts = glob_wildcards("Reads/{read}.{ext}");
threads = 8;

ML_PRESETS = {
	"min5":  "-k 5 -w 11",
	"min7": "-k 7 -w 13",
}

rule run:
	input:
		expand(expand("Results/{read}.{ext}|minimap2|{{tp}}|0|1|2|REPORT/",
					  zip, read=reads, ext=exts),
			   tp = ["ONT", "PB"]),
		expand(expand("Results/{read}.{ext}|miniLyndon|{{preset}}|CFL_ICFL|{{comb}}|{{segment_size}}|REPORT/",
					  zip, read=reads, ext=exts),
			   segment_size=[-1,100,300,500],
			   preset=ML_PRESETS.keys(),
			   comb=["NONCOMB", "COMB"]),
		expand(expand("Results/{read}.{ext}|miniLyndon|{{preset}}|CFL_ICFL_R|{{comb}}|{{recursive_size}}|REPORT/",
					  zip, read=reads, ext=exts),
			   recursive_size=[10,25,50],
			   preset = ML_PRESETS.keys(),
			   comb = ["NONCOMB", "COMB"])

rule miniLyndon:
	input:
		"Reads/{read}.{ext}"
	output:
		"Results/{read}.{ext}|miniLyndon|{preset}|{factorization}|{comb}|{size}.paf"
	params:
		segment_recursive_size = lambda wildcards: "-s "+wildcards.size if wildcards.factorization == "CFL_ICFL" else "-r "+wildcards.size,
		comb_value=lambda wildcards: 1 if wildcards.comb == "COMB" else 0,
		preset_params= lambda wildcards: ML_PRESETS[wildcards.preset],
		reference=reference,
		threads=threads
	threads: threads
	shell:
		"""
		echo "1---------------------------------MiniLyndon---------------------------------1\n" &&
		{{ time ../miniLyndon/bin/fingerprint -f "{wildcards.factorization}" -p "Reads/" -a "{wildcards.read}.{wildcards.ext}" -n {threads} {params.segment_recursive_size} -c {params.comb_value} | \
		../miniLyndon/bin/minimizer_demo -t {threads} {params.preset_params} | \
		../miniLyndon/bin/postprocessing "{input}" > "{output}" ; }} 2>&1 | tee -a "Results/{wildcards.read}.{wildcards.ext}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}-benchmark.txt"
		"""

rule minimap:
	input:
		"Reads/{read}.{ext}"
	output:
		"Results/{read}.{ext}|minimap2|{tp}|{factorization}|{comb}|{size}.paf"
	params:
		read_type = lambda wildcards: "-x ava-ont" if wildcards.tp == "ONT" else "-x ava-pb",
		reference=reference
	threads: threads
	shell:
		"""
		echo "1---------------------------------Minimap2---------------------------------1\n" &&
		{{ time ./minimap2/minimap2 {params.read_type} -t {threads} "{input}" "{input}" > "{output}" ; }} 2>&1 | \
		tee -a "Results/{wildcards.read}.{wildcards.ext}|minimap2|{wildcards.tp}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}-benchmark.txt"
		"""

rule miniasm:
	input:
		read="Reads/{read}.{ext}",
		paf="Results/{read}.{ext}|{tool}|{preset}|{factorization}|{comb}|{size}.paf"
	output:
		gfa="miniasm/{read}.{ext}|{tool}|{preset}|{factorization}|{comb}|{size}.gfa",
		fa="miniasm/{read}.{ext}|{tool}|{preset}|{factorization}|{comb}|{size}.fa"
	shell:
		"""
		echo "2---------------------------------Miniasm---------------------------------2\n" &&
		./miniasm/miniasm -f "{input.read}" "{input.paf}" > "{output.gfa}" &&
		awk '/^S/{{print \">\" $2 \"\\n\" $3}}' {output.gfa} | fold > {output.fa}
		"""

rule quast:
	input:
		"miniasm/{read}.{ext}|{tool}|{preset}|{factorization}|{comb}|{size}.fa"
	output:
		directory("Results/{read}.{ext}|{tool}|{preset}|{factorization}|{comb}|{size}|REPORT")
	params:
		reference=reference
	threads: threads
	shell:
		"""
		echo "3---------------------------------Quast---------------------------------3\n" &&
		./quast/quast.py --threads {threads} "{input}" -r {params.reference} -o "{output}" &&
		echo "4---------------------------------Cleanup---------------------------------4\n" &&
		mv -v -f "Results/{wildcards.read}.{wildcards.ext}|{wildcards.tool}|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}-benchmark.txt" "./Results/{wildcards.read}.{wildcards.ext}|{wildcards.tool}|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}|REPORT/benchmark.txt"
		"""

rule clean:
	shell:
		'''
		rm -f -r Results
		rm -f ./miniasm/*.gfa
		rm -f ./miniasm/*.fa
		'''

# snakemake --config reference="/home/omarm/Desktop/Paper/Sperimentazione_Lyndon/Sperimentazione/Ecoli_K12_DH10B.fasta" -s make.smk -f run --cores 16 --keep-going
# snakemake --config reference="/home/omarm/Desktop/Paper/Sperimentazione_Lyndon/Sperimentazione/Ecoli_K12_DH10B.fasta" -s make.smk -f clean --cores 1
