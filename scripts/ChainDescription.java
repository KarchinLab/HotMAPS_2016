import java.util.*;
import java.io.*;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.StructureException;

// Puts a String array and a AttributeField object into one single object for comparator
// Perhaps not the best you of a class but it avoids type casting issues when inheriting from comparator
class LineData {

    ArrayList<String> line;

    public LineData(String[] splitted_line) {
        this.line = new ArrayList( Arrays.asList( splitted_line ) );
    }

    // getter for string array
    public ArrayList<String> getLine(){
        return this.line;
    }
}

// Read in tab delimited gtf file
class TabDelimReader {

    ArrayList<LineData> lines = new ArrayList<LineData>();
    public TabDelimReader(String fname){
        try {
            String[] tmp_line;
            Scanner scan = new Scanner(new File(fname));
            while(scan.hasNextLine()){
                tmp_line = scan.nextLine().split("\t");  // hold line data temporarily
                this.lines.add(new LineData(tmp_line));
            }
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    public ArrayList<LineData> getLines(){
        return this.lines;
    }
}


class SortOutputComparator implements Comparator<Object> {

    public int compare(Object o1, Object o2){
        LineData l1 = (LineData) o1;
        LineData l2 = (LineData) o2;
        ArrayList<String> line1 = l1.getLine();
        ArrayList<String> line2 = l2.getLine();

        // compare chromosome names
        String pdbId1 = (String) line1.get(0);
        int pdbIdCompare = pdbId1.compareTo((String) line2.get(0));
        String description1 = (String) line1.get(3);
        int descriptCompare = pdbId1.compareTo((String) line2.get(3));
        boolean headerCompare1 = line1.get(0).equals(new String("PDBId"));
        boolean headerCompare2 = line2.get(0).equals(new String("PDBId"));

        // compare gene/tx ids
        try{
            // comparison logic
            if (headerCompare1){
                return -1; 
            } else if (headerCompare2){
                return 1; 
            } else if(pdbIdCompare < 0){
                return -1;
            } else if (pdbIdCompare > 0){
                return 1;
            } else if (descriptCompare < 0) {
                return -1;
            } else if (descriptCompare > 0) {
                return 1;
            } else {
                return 0;
            }
        } catch(Exception e){
            e.printStackTrace();
            System.exit(1);
            return 0;
        }
    }
}


// SortGtf is the main class that drives sorting of a GTF by [seqname, gene id, tx id, start, end]
public class ChainDescription {

    // The join method acts just like python's string.join()
    public static String join(String[] fields, String separator){
        StringBuilder stringBuilder = new StringBuilder();
        for(int i = 0; i < fields.length; i++){
            stringBuilder.append(fields[i]);
            if(i < fields.length - 1){
                stringBuilder.append(separator);
            }
        }
        return stringBuilder.toString();
    }


    public static void main(String[] args){
        // complain about args
        if (args.length != 2){
            System.out.println("ChainDescription.java requires exactly two arguments.");
            System.out.println("Usage: java ChainDescription input.txt output.txt");
            System.out.println("ChainDescription adds the chain description to the PDB info file");
            System.exit(1);
        }

        // Read in the PDB info file and then sort
        TabDelimReader inputReader = new TabDelimReader(args[0]);
        ArrayList<LineData> lines = inputReader.getLines();  // read in file (tab-separated)
	// list for valid line
	ArrayList<LineData> found_lines = new ArrayList<LineData>();
	// separate list for missing files
	ArrayList<LineData> missing_lines = new ArrayList<LineData>();
        // add header column
        lines.get(0).getLine().add("ChainDescription");
        LineData header = lines.get(0);

        // iterate through all lines
        int numLines = lines.size(), i = 0;
        LineData myLine = null;
        String pdbFilePath = null;

        for (i=1;i<numLines;i++){
            myLine = lines.get(i);

            // check if it is case of missing pdb file                                                                         
	    if (myLine.getLine().size() < 4) {
		//System.out.println("Here");
		missing_lines.add(myLine);
	    }
            
	    // pdb file exists 
	    else if (!myLine.getLine().get(3).equals("")) {
                pdbFilePath = myLine.getLine().get(3);

                // Read in PDB file
                try {
                    PDBFileReader pdbReader = new PDBFileReader();
                    Structure structure = null;
                    structure = pdbReader.getStructure(pdbFilePath);
                    Chain chain = null;
                    
                    // get specific chain
                    try{
                        chain = structure.findChain(myLine.getLine().get(1));
                    } catch(StructureException s1){
                        System.out.println(myLine);
                        s1.printStackTrace();
                        return;
                    }

                    String chainDesc = chain.toString().split("\n")[1].replace(" ", "_").trim(); 
                    myLine.line.add(chainDesc);
		    found_lines.add(myLine);
                    //}
                } catch (IOException e1) {
                    e1.printStackTrace();
                    System.out.println("PDB File Path:");
                    System.out.println(pdbFilePath);
                    return;
                } 
	    }
	    // otherwise it is case of homology model
	    else {
                myLine.line.add("");
		found_lines.add(myLine);
            }
        }


        // sort the output file
        Object[] data = found_lines.toArray();
        Arrays.sort(data, new SortOutputComparator()); // sort the data
	
	// missing data handling
	Object[] missing_data = missing_lines.toArray();

        // output the new pdb info file
        try{
            Formatter output=new Formatter(args[1]);
	    	    
            String[] dataLine = new String[6];
            LineData tmp;
	    
	    // add the header
	    tmp = (LineData) header;
	    dataLine = tmp.getLine().toArray(dataLine);
	    output.format("%s\n", ChainDescription.join(dataLine, "\t"));
	    
	    // add the found pdb_ids
	    for(i=0; i<data.length;i++){
                tmp = (LineData) data[i];
                dataLine = tmp.getLine().toArray(dataLine);
                output.format("%s\n", ChainDescription.join(dataLine, "\t"));
            }

	    // add the missing pdb_ids
	    dataLine = new String[6];
	    for (i=0; i<missing_data.length; i++) {
		tmp = (LineData) missing_data[i];
		dataLine = tmp.getLine().toArray(dataLine);
		dataLine[3] = "";
		dataLine[4] = "";
		dataLine[5] = "";
		output.format("%s\n", ChainDescription.join(dataLine, "\t"));
	    }
            output.close();
        } catch(FileNotFoundException e){
            e.printStackTrace();
        }

    }
}
