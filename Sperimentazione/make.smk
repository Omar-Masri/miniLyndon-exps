reference = config["reference"]

reads, = glob_wildcards("Reads/{read}.fa");

threads = 8;

rule run:
	input:
		expand("Results/{read}|miniLyndon|{preset}|CFL_ICFL|{comb}|{segment_size}.paf",
			   read=reads,
			   segment_size=[-1,100,300,500],
			   preset=["min5", "min7"],
			   comb=["NONCOMB", "COMB"]),
		expand("Results/{read}|miniLyndon|{preset}|CFL_ICFL_R|{comb}|{recursive_size}.paf",
			   read=reads,
			   recursive_size=[10,25,50],
			   preset = ["min5", "min7"],
			   comb = ["NONCOMB", "COMB"])

rule CFL_ICFL:
	input:
		read="Reads/{read}.fa"
	output:
		"Results/{read}|miniLyndon|{preset}|{factorization}|{comb}|{size}.paf"
	params:
		segment_recursive_size = lambda wildcards: "-s "+wildcards.size if wildcards.factorization == "CFL_ICFL" else "-r "+wildcards.size,
		comb_value=lambda wildcards: 1 if wildcards.comb == "COMB" else 0,
		preset_params=lambda wildcards: "-k 7 -w 13" if wildcards.preset == "min7" else "-k 5 -w 11",
		reference=reference,
		threads=threads
	shell:
		'''
		echo "1---------------------MiniLyndon---------------------1\n" &&
		{{ time ../miniLyndon/bin/fingerprint -f "{wildcards.factorization}" -p "Reads/" -a "{wildcards.read}.fa" -n {params.threads} {params.segment_recursive_size} -c {params.comb_value} | \
		../miniLyndon/bin/minimizer_demo -t {params.threads} {params.preset_params} | \
		../miniLyndon/bin/postprocessing "Reads/{wildcards.read}.fa" > "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}.paf" ; }} 2>&1 | tee -a "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}-benchmark.txt" &&
		echo "2---------------------Miniasm---------------------2\n" &&
		./miniasm/miniasm -f Reads/{wildcards.read}.fa "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}.paf" > "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}.gfa" &&
		awk '/^S/{{print \">\" $2 \"\\n\" $3}}' "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}.gfa" | fold > "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}.fa" &&
		echo "3---------------------Quast---------------------3\n" &&
		./quast/quast.py --threads {params.threads} "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}.fa" -r {params.reference} -o ./quast/quast_test_output &&
		echo "4---------------------Cleanup---------------------4\n" &&
		mv -v -f ./quast/quast_test_output/ "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}|REPORT/" &&
		mv -v -f  "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}-benchmark.txt" "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|{wildcards.factorization}|{wildcards.comb}|{wildcards.size}|REPORT/benchmark.txt" || echo "FAILS"
		'''

rule clean:
	shell:
		'''
		rm -f -r Results
		rm -f ./miniasm/*.gfa
		rm -f ./miniasm/*.fa
		'''

# snakemake --config reference="/home/omarm/Desktop/Paper/Sperimentazione_Lyndon/Sperimentazione/Ecoli_K12_DH10B.fasta" -s make.smk -f run --cores 1
# snakemake --config reference="/home/omarm/Desktop/Paper/Sperimentazione_Lyndon/Sperimentazione/Ecoli_K12_DH10B.fasta" -s make.smk -f clean --cores 1
