/*
 * Code for correcting SV sequence, particularly INS and INV
 * basically taking the code from PhaseSVs Java, but only keep the stuff for INV population
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
	static String USAGE = "splicephase.pl phased.vcf "
    + "sniffles.vcf loadreads.hairs spliced.vcf ref.fa\n";

    public static void main(String[] args) throws Exception
	{
		if(args.length < 5)
		{
			System.out.println(USAGE);
			System.exit(1);
		}
		else
		{
			// Check that the input files all exist
			
			if(!new File(args[0]).exists())
			{
				System.out.println("Invalid phased SNP file: " + args[0]);
				System.exit(1);
			}
			
			if(!new File(args[1]).exists())
			{
				System.out.println("Invalid Sniffles call file: " + args[1]);
				System.exit(1);
			}
			
			if(!new File(args[2]).exists())
			{
				System.out.println("Invalid HapCUT2 fragment file: " + args[1]);
				System.exit(1);
			}
			
			if(!new File(args[4]).exists())
			{
				System.out.println("Invalid reference file: " + args[4]);
				System.exit(1);
			}
			
			runSplicing(args[0], args[1], args[2], args[3], args[4]);
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
		Process child = Runtime.getRuntime().exec(faidxCommand);
		InputStream seqStream = child.getInputStream();
		Scanner seqInput = new Scanner(seqStream);
		
		// Make sure it produced an actual output
		if(!seqInput.hasNext())
		{
			seqInput.close();
			throw new Exception("samtools faidx did not produce an output: " + faidxCommand);
		}
		// Read in and ignore sequence name
		seqInput.next();
		
		// Make sure there's a sequence
		if(!seqInput.hasNext())
		{
			seqInput.close();
			throw new Exception("samtools faidx produced a sequence name but not an actual sequence: " + faidxCommand);
		}
		
		// Concatenate all lines of the output sequence
		StringBuilder res = new StringBuilder("");
		while(seqInput.hasNext())
		{
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

	
}