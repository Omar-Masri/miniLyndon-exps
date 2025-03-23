genomes, g_exts = glob_wildcards("Reference_Genomes/{read}.{ext}");

reads, r_exts = glob_wildcards("Reads/{read}.{ext}");
threads = 8;

ML_PRESETS = {
	"min5":  "-k 5 -w 11",
	"min7": "-k 7 -w 13",
}

rule run:
	input:
		expand(
			expand(
				expand(
					"Results/{{reference}}.{{g_ext}}/{read}.{r_ext}|minimap2|{{{{tp}}}}|0|1|2|REPORT/",
					zip, read=reads, r_ext=r_exts),
				zip, reference=genomes, g_ext=g_exts),
			tp = ["ONT", "PB"]),
		expand(
			expand(
				expand(
					"Results/{{reference}}.{{g_ext}}/{read}.{r_ext}|miniLyndon|{{{{preset}}}}|CFL_ICFL_R|{{{{comb}}}}|{{{{segment_size}}}}|REPORT/",
					zip, read=reads, r_ext=r_exts),
				zip, reference=genomes, g_ext=g_exts),
			segment_size=[-1,100,300,500],
			preset = ML_PRESETS.keys(),
			comb = ["NONCOMB", "COMB"]),
		expand(
			expand(
				expand(
					"Results/{{reference}}.{{g_ext}}/{read}.{r_ext}|miniLyndon|{{{{preset}}}}|CFL_ICFL_R|{{{{comb}}}}|{{{{recursive_size}}}}|REPORT/",
					zip, read=reads, r_ext=r_exts),
				zip, reference=genomes, g_ext=g_exts),
			recursive_size=[10,25,50],
			preset = ML_PRESETS.keys(),
			comb = ["NONCOMB", "COMB"])

rule miniLyndon:
	input:
		"Reads/{read}.{r_ext}"
	output:
		paf="Results/{reference}.{g_ext}/{read}.{r_ext}|miniLyndon|{preset}|{factorization}|{comb}|{size}.paf",
		txt="Results/{reference}.{g_ext}/{read}.{r_ext}|miniLyndon|{preset}|{factorization}|{comb}|{size}-benchmark.txt"
	params:
		segment_recursive_size = lambda wildcards: "-s "+wildcards.size if wildcards.factorization == "CFL_ICFL" else "-r "+wildcards.size,
		comb_value=lambda wildcards: 1 if wildcards.comb == "COMB" else 0,
		preset_params= lambda wildcards: ML_PRESETS[wildcards.preset]
	threads: threads
	shell:
		"""
		echo "1---------------------------------MiniLyndon---------------------------------1\n" &&
		{{ time ../miniLyndon/bin/fingerprint -f "{wildcards.factorization}" -p "Reads/" -a "{wildcards.read}.{wildcards.r_ext}" -n {threads} {params.segment_recursive_size} -c {params.comb_value} | \
		../miniLyndon/bin/minimizer_demo -t {threads} {params.preset_params} | \
		../miniLyndon/bin/postprocessing "{input}" > "{output.paf}" ; }} 2>&1 | tee -a "{output.txt}"
		"""

rule minimap:
	input:
		"Reads/{read}.{r_ext}"
	output:
		paf="Results/{reference}.{g_ext}/{read}.{r_ext}|minimap2|{tp}|{factorization}|{comb}|{size}.paf",
		txt="Results/{reference}.{g_ext}/{read}.{r_ext}|minimap2|{tp}|{factorization}|{comb}|{size}-benchmark.txt"
	params:
		read_type = lambda wildcards: "-x ava-ont" if wildcards.tp == "ONT" else "-x ava-pb",
	threads: threads
	shell:
		"""
		echo "1---------------------------------Minimap2---------------------------------1\n" &&
		{{ time minimap2 {params.read_type} -t {threads} "{input}" "{input}" > "{output.paf}" ; }} 2>&1 | \
		tee -a "{output.txt}"
		"""

rule miniasm:
	input:
		read="Reads/{read}.{r_ext}",
		paf="Results/{reference}.{g_ext}/{read}.{r_ext}|{tool}|{preset}|{factorization}|{comb}|{size}.paf"
	output:
		gfa="miniasm/{reference}.{g_ext}/{read}.{r_ext}|{tool}|{preset}|{factorization}|{comb}|{size}.gfa",
		fa="miniasm/{reference}.{g_ext}/{read}.{r_ext}|{tool}|{preset}|{factorization}|{comb}|{size}.fa"
	threads: 1
	shell:
		"""
		echo "2---------------------------------Miniasm---------------------------------2\n" &&
		miniasm -f "{input.read}" "{input.paf}" > "{output.gfa}" &&
		awk '/^S/{{print \">\" $2 \"\\n\" $3}}' "{output.gfa}" | fold > "{output.fa}"
		"""

rule quast:
	input:
		genome="Reference_Genomes/{reference}.{g_ext}",
		fa="miniasm/{reference}.{g_ext}/{read}.{r_ext}|{tool}|{preset}|{factorization}|{comb}|{size}.fa"
	output:
		directory=directory("Results/{reference}.{g_ext}/{read}.{r_ext}|{tool}|{preset}|{factorization}|{comb}|{size}|REPORT"),
		benchmark="Results/{reference}.{g_ext}/{read}.{r_ext}|{tool}|{preset}|{factorization}|{comb}|{size}|REPORT/benchmark.txt"
	threads: 4
	shell:
		"""
		echo "3---------------------------------Quast---------------------------------3\n" &&
		quast --threads {threads} "{input.fa}" -r "{input.genome}" -o "{output.directory}" &&
		echo "4---------------------------------Cleanup---------------------------------4\n" &&
		mv -v -f "Results/{wildcards.reference}.{wildcards.g_ext}/{wildcards.read}.{wildcards.r_ext}|{wildcards.tool}|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}-benchmark.txt" "{output.benchmark}"
		"""

rule clean:
	shell:
		'''
		rm -f -r Results
		rm -f -r ./miniasm/*
		'''

# snakemake -s make.smk -f run --cores 2 --keep-going
# snakemake -s make.smk -f clean --cores 1 --keep-going

# TODO BENCHMARK:
