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
		"Results/{read}|miniLyndon|{preset}|CFL_ICFL|{comb}|{segment_size}.paf"
	params:
		comb_value= lambda wildcards: 1 if wildcards.comb == "COMB" else 0,
		preset_params=lambda wildcards: "-k 7 -w 13" if wildcards.preset == "min7" else "-k 5 -w 11",
		reference=reference
	threads: threads
	shell:
		'''
		../miniLyndon/bin/fingerprint -f "CFL_ICFL" -p "Reads/" -a "{wildcards.read}.fa" -n {threads} -s {wildcards.segment_size} -c {params.comb_value} | ../miniLyndon/bin/minimizer_demo -t {threads} {params.preset_params} | ../miniLyndon/bin/postprocessing "Reads/{wildcards.read}.fa" > "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}.paf" &&
		./miniasm/miniasm -f Reads/{wildcards.read}.fa "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}.paf" > "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}.gfa" &&
		awk '/^S/{{print \">\" $2 \"\\n\" $3}}' "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}.gfa" | fold > "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}.fa" &&
		./quast/quast.py "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}.fa" -r {params.reference} -o quast_test_output &&
		mv ./quast/quast_test_output/report.pdf "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}|REPORT.pdf" &&
		mv ./quast/quast_test_output/report.txt "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL|{wildcards.comb}|{wildcards.segment_size}|REPORT.txt" || echo "FAILS"
		'''

rule CFL_ICFL_R:
	input:
		read="Reads/{read}.fa"
	output:
		"Results/{read}|miniLyndon|{preset}|CFL_ICFL_R|{comb}|{recursive_size}.paf"
	params:
		comb_value= lambda wildcards: 1 if wildcards.comb == "COMB" else 0,
		preset_params=lambda wildcards: "-k 7 -w 13" if wildcards.preset == "min7" else "-k 5 -w 11",
		reference=reference
	threads: threads
	shell:
		'''
		../miniLyndon/bin/fingerprint -f "CFL_ICFL_R" -p "Reads/" -a "{wildcards.read}.fa" -n {threads} -r {wildcards.recursive_size} -c {params.comb_value} | ../miniLyndon/bin/minimizer_demo -t {threads} {params.preset_params} | ../miniLyndon/bin/postprocessing "Reads/{wildcards.read}.fa" > "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}.paf" &&
		./miniasm/miniasm -f Reads/{wildcards.read}.fa "Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}.paf" > "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}.gfa" &&
		awk '/^S/{{print \">\" $2 \"\\n\" $3}}' "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}.gfa" | fold > "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}.fa" &&
		./quast/quast.py "./miniasm/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}.fa" -r {params.reference} -o quast_test_output &&
		mv ./quast/quast_test_output/report.pdf "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}|REPORT.pdf" &&
		mv ./quast/quast_test_output/report.txt "./Results/{wildcards.read}|miniLyndon|{wildcards.preset}|CFL_ICFL_R|{wildcards.comb}|{wildcards.recursive_size}|REPORT.txt" || echo "FAILS"
		'''

rule clean:
	shell:
		'''
		rm -f -r Results
		rm -f ./miniasm/*.gfa
		rm -f ./miniasm/*.fa
		'''
