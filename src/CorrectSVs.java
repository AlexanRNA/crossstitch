/*
 * Code for correcting SV sequence, particularly INS and INV
 * basically taking the code from CorrectSVs Java, but only keep the stuff for INV population
 * Ignore any output from the report, since it will be wrong! It is the part of the code I left 
 * unchanged 
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;


public class CorrectSVs {
	// Maximum length of an SV to include
	static int MAX_SV_LEN = 100000;

	// The path to the samtools executable
	static String SAMTOOLS_PATH = "samtools";

	// Usage message if running the program incorrectly
	static String USAGE = " phased.vcf "
			+ "sniffles.vcf spliced.vcf ref.fa\n";

	public static void main(String[] args) throws Exception
	{
		if(args.length < 3)
		{
			System.out.println(USAGE);
			System.exit(1);
		}
		else
		{ // check if the file exist
			if(!new File(args[0]).exists())
			{
				System.out.println("Invalid Sniffles call file: " + args[0]);
				System.exit(1);
			}

			if(!new File(args[1]).exists())
			{
				System.out.println("Invalid HapCUT2 fragment file: " + args[1]);
				System.exit(1);
			}

			if(!new File(args[2]).exists())
			{
				System.out.println("Invalid reference file: " + args[2]);
				System.exit(1);
			}

			runCorrection(args[0], args[1], args[2]);
		}
	}

	/*
	 * Reverse complement of a String
	 */
	static char[] froms = new char[] {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
	static char[] tos = new char[] {'T', 'G', 'C', 'A', 't', 'g', 'c', 'a'};
	static String revComp(String s)
	{
		int n = s.length();
		char[] res = new char[n];
		for(int i = 0; i<n; i++)
		{
			char c = s.charAt(i);
			char to = c;
			for(int j = 0; j<froms.length; j++)
			{
				if(c == froms[j])
				{
					to = tos[j];
					break;
				}
			}
			res[n-1-i] = to;
		}
		return new String(res);
	}

	/*
	 * Queries a genomic substring - runs samtools faidx <genomeFile> chr:startPos-endPos
	 */
	static String genomeSubstring(String refFn, String chr, long startPos, long endPos) throws Exception
	{
		if(startPos > endPos)
		{
			return "";
		}
		String faidxCommand = String.format("%s faidx %s %s:%d-%d", SAMTOOLS_PATH, refFn, chr, startPos, endPos);
		ProcessBuilder processBuilder = new ProcessBuilder(faidxCommand.split(" "));
		Process process = processBuilder.start();

		InputStream seqStream = process.getInputStream();
		Scanner seqInput = new Scanner(seqStream);

		// Make sure it produced an actual output
		if (!seqInput.hasNext()) {
			seqInput.close();
			throw new Exception("samtools faidx did not produce an output: " + faidxCommand);
		}
		// Read in and ignore sequence name
		seqInput.next();

		// Make sure there's a sequence
		if (!seqInput.hasNext()) {
			seqInput.close();
			throw new Exception("samtools faidx produced a sequence name but not an actual sequence: " + faidxCommand);
		}

		// Concatenate all lines of the output sequence
		StringBuilder res = new StringBuilder("");
		while (seqInput.hasNext()) {
			res.append(seqInput.next());
		}
		seqInput.close();

		return res.toString();
	}

	/*
	 * Gets the reference sequence for a variant using samtools
	 */
	static String getVariantSeq(String chr, long startPos, long svLen, String refFn) throws Exception
	{
		return genomeSubstring(refFn, chr, startPos, startPos + svLen);
	}
	static void runCorrection(//String snpsFn,
							  String svsFn, String vcfOfn,
							  String refFn) throws Exception {
		// Generate filenames for phasing information files based on output VCF
		String readPhaseOfn = vcfOfn + ".readphase";
		String svPhaseOfn = vcfOfn + ".svphase";
		String svPhaseDetailsOfn = vcfOfn + ".svphase.details";

		// Initialize PrintWriters for output files
		PrintWriter out = new PrintWriter(new File(vcfOfn));
		PrintWriter readPhaseOut = new PrintWriter(new File(readPhaseOfn));
		PrintWriter svPhaseOut = new PrintWriter(new File(svPhaseOfn));
		PrintWriter svPhaseDetailsOut = new PrintWriter(new File(svPhaseDetailsOfn));

		// Initialize Scanners for input files
		//Scanner snpsInput = new Scanner(new FileInputStream(new File(snpsFn)));
		Scanner svsInput = new Scanner(new FileInputStream(new File(svsFn)));
		// Scanner hairsInput = new Scanner(new FileInputStream(new File(hairsFn)));

		// Now load the sniffles calls
		System.err.println("Loading Sniffles calls");
		ArrayList<String> svHeader = new ArrayList<String>();
		int svCount = 0;

		// Information to report about each type
		HashMap<String, CorrectSVs.Report> typeReports = new HashMap<String, CorrectSVs.Report>();
		HashMap<String, CorrectSVs.Read> readsToPhase = new HashMap<String, CorrectSVs.Read>();
		TreeMap<String, TreeMap<Integer, CorrectSVs.Variant>> svData = new TreeMap<String, TreeMap<Integer, CorrectSVs.Variant>>();

		// Read SVs one at a time
		while (svsInput.hasNext()) {

			String line = svsInput.nextLine();
			if (line.startsWith("#")) {
				svHeader.add(line);
				 // continue;
			} else {

				svCount++;

				CorrectSVs.Variant v = new CorrectSVs.Variant(line);
				v.fillSvFields();

				// Update counter for the type
				if (!typeReports.containsKey(v.type)) {
					typeReports.put(v.type, new CorrectSVs.Report());
				}
				typeReports.get(v.type).total++;

				// Update counters for each read supporting this SV
				for (String s : v.readNames) {
					if (!readsToPhase.containsKey(s)) {
						readsToPhase.put(s, new CorrectSVs.Read());
					}
					readsToPhase.get(s).count++;
				}

				// Update map of (chromosome, position) pair to SV
				if (!svData.containsKey(v.chromosome)) {
					svData.put(v.chromosome, new TreeMap<Integer, CorrectSVs.Variant>());
				}
				svData.get(v.chromosome).put(v.pos, v);
			}
		}

				svsInput.close();

			System.err.printf("Loaded Sniffles calls: %d variants involving %d reads\n", svCount, readsToPhase.size());
			for (String s : typeReports.keySet()) {
				System.err.println(s + ": " + typeReports.get(s).total);
			}

			// Log some information about read phasing
			//readPhaseOut.println("#READID\tNUMSV\t|\tHAP1\tHAP2\t| HAP HAPR | SNPS\n");
			ArrayList<String> sortedReadNames = new ArrayList<String>();
			sortedReadNames.addAll(readsToPhase.keySet());
			Collections.sort(sortedReadNames);
			for (String s : sortedReadNames) {
				CorrectSVs.Read r = readsToPhase.get(s);
				//double hap1Prop = (r.hap1 + r.hap2 > 0) ? 100.0 * r.hap1 / (r.hap1 + r.hap2) : 0;
				//String hap = r.hap1 >= r.hap2 ? "hapA" : "hapB";
				//readPhaseOut.printf("%s\t%d\t|\t%d\t%d\t| %s %7.02f  |%s\n",
				//		s, r.count, r.hap1, r.hap2, hap, hap1Prop, r.snps.toString());
			}
			readPhaseOut.close();

			//System.err.printf("Scanned fragments: %d lines with %d involved in SVs\n", hairsLines, svHairs);

			// Process Sniffles SVs and phase them based on their reads
			svPhaseOut.println("chr:pos:genotype\ttype\tsvlen\tseqlen\t|\tnumreads\t|\tincludesv");
			svPhaseDetailsOut.println("chr:pos:genotype\ttype\tsvlen\tseqlen\t|\tincludesv");
			int snifflesCount = 0, reportedSvs = 0, unphasedSvs = 0, svLenErr = 0, genotypeErr = 0;

			// Iterate over the SVs one chromosome at a time
			for (String chr : svData.keySet()) {
				TreeMap<Integer, CorrectSVs.Variant> svChrData = svData.get(chr);

				// Loop over the positions with variants on this chromosome
				for (int pos : svChrData.keySet()) {
					snifflesCount++;
					CorrectSVs.Variant v = svChrData.get(pos);


					String oldGenotype = v.genotype;

					if (v.type.equals("INS") || v.type.equals("DEL") || v.type.equals("INV") || v.type.equals("DUP")) {
						if (Math.abs(v.svlen) <= MAX_SV_LEN) {

							// This is the sequence length of the SV and is used for logging
							int seqLength = 0;
							if (v.seq != null) {
								seqLength = v.seq.length();
							}

							for (String readName : v.readNames) {
								svPhaseDetailsOut.println("== " + readName);
							}

							// Whether or not to include this SV - can be set to false based on certain criteria
							boolean includeSv = true;

							if (v.type.equals("DEL")) {
								continue;
							} else if (v.type.equals("DUP")) {

								//StringBuilder refBuilder = new StringBuilder("");
								if (v.alt.equals("<DUP>")){

									v.alt = getVariantSeq(chr, (long) pos, (long) v.svlen, refFn);
									// v.alt = refBuilder.toString();
								}
							}else if (v.type.equals("INS")) {
								StringBuilder refBuilder = new StringBuilder("");
								for (int i = 0; i < v.svlen; i++) {
									refBuilder.append("N");
								}
								if (v.alt.equals("<INS>")){
									// if the alt seq is <INS>, replace it with Ns where the number of Ns corresponds
									//to the length of INS
									v.alt = refBuilder.toString();
								}
							} else if (v.type.equals("INV")) {
								StringBuilder refBuilder = new StringBuilder("");
								for (int i = 0; i < v.svlen; i++) {
									refBuilder.append("X");
								}
								v.ref = refBuilder.toString();
								v.alt = revComp(getVariantSeq(chr, (long) pos, (long) v.svlen, refFn));
							}

							// Log this SV to both logging files
							svPhaseOut.printf("%s:%d:%s\t<%s>\t%d\t%d\t| \t%s\t%d\n",
									chr, pos, oldGenotype, v.type, v.svlen, seqLength, v.readNames.length, includeSv ? 1 : 0);

							svPhaseDetailsOut.printf("%s:%d:%s\t<%s>\t%d\t%d\t| \t%s\t%d\n",
									chr, pos, oldGenotype, v.type, v.svlen, seqLength, v.readNames.length,
									includeSv ? 1 : 0);

							// If we include this in our phasing, mark it as phased and update the genotype
							if (includeSv) {

								v.sample = oldGenotype + v.sample.substring(3);
								svData.get(chr).put(pos, v);
								reportedSvs++;

								typeReports.get(v.type).reported++;


								typeReports.get(v.type).unphased++;
								unphasedSvs++;

							}

						} else {
							System.err.println("ERROR: extreme SV length reported: " + v.type + " " + v.svlen);
							svLenErr++;
						}
					}
				}
			}

			svPhaseOut.close();
			svPhaseDetailsOut.close();

			// Output some SV phasing statistics
			System.err.printf("Reported %d, unphased %d of %d attempted.  svlenerr: %d, genotypeerr: %d\n",
					reportedSvs, unphasedSvs, snifflesCount, svLenErr, genotypeErr);
			System.err.println("type all reported phased unphased:");

			for (String s : typeReports.keySet()) {
				CorrectSVs.Report r = typeReports.get(s);
				System.err.println(s + " " + r.total + " " + r.reported + " " + r.unphased);
			}

			System.err.println();

			// Print header to output file
			for (String headerLine : svHeader) {
				out.println(headerLine);
			}

			// Print final variants to output file
			for (String chr : svData.keySet()) {
				for (int pos : svData.get(chr).keySet()) {
					CorrectSVs.Variant v = svData.get(chr).get(pos);
					out.println(v);
				}
			}

			out.close();

	}

	/*
	 * Information for each individual read
	 */
	static class Read
	{
		// The number of SVs it's included in
		int count;


		StringBuilder snps;
		Read()
		{
			count = 0;
			//	hap1 = 0;
			//	hap2 = 0;
			snps = new StringBuilder("");
		}
	}

	/*
	 * Information to be reported for each type of variant
	 */
	static class Report
	{
		int total;
		int reported;
		//int phased;
		int unphased;
		Report()
		{
			total = reported  = unphased = 0;
		}
	}

	/*
	 * Information for a variant
	 */
	static class Variant
	{
		String chromosome;
		int pos;
		String id, ref, alt;
		String qual, filter, info;
		String format, sample;
		String genotype;
		String[] readNames;
		String type;
		String seq;
		int svlen;
		int svend;
		Variant(String line)
		{
			String[] tokens = line.split("\t");
			chromosome = tokens[0];
			pos = Integer.parseInt(tokens[1]);
			id = tokens[2];
			ref = tokens[3];
			alt = tokens[4];
			qual = tokens[5];
			filter = tokens[6];
			info = tokens[7];
			format = tokens[8];
			sample = tokens.length > 9 ? tokens[9] : "";
			if(sample.length() > 0)
			{
				genotype = sample.substring(0, sample.indexOf(":"));
			}
		}

		/*
		 * Populates list of supporting reads from the RNAMES INFO field
		 */
		void fillSvFields()
		{
			String[] tokens = info.split(";");

			for(String token : tokens)
			{
				int equalsIdx = token.indexOf('=');
				if(equalsIdx != -1)
				{
					String key = token.substring(0, equalsIdx);
					String val = token.substring(1 + equalsIdx);
					if(key.equals("RNAMES"))
					{
						readNames = val.split(",");
					}
					else if(key.equals("SEQ"))
					{
						seq = val;
					}
					else if(key.equals("SVLEN"))
					{
						svlen = Math.abs(Integer.parseInt(val));
					}
					else if(key.equals("END"))
					{
						svend = Math.abs(Integer.parseInt(val));
					}
					else if(key.equals("SVTYPE"))
					{
						type = val;
					}
				}
			}
			if(type == null && alt.contains("<"))
			{
				type = alt.substring(1, alt.length()-1);
			}
		}

		/*
		 * VCF formatted line (with no new line at the end)
		 */
		public String toString()
		{
			return chromosome + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual
					+ "\t" + filter + "\t" + info + "\t" + format + (sample.length() == 0 ? "" : ("\t" + sample));
		}
	}

}