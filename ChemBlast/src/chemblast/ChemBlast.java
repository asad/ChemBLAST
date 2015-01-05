/*
 // This Java source file is copyright (C) 2014 by Syed Asad Rahman. All rights
 // reserved. For further information, contact the author, Syed Asad Rahman, at
 // asad@ebi.ac.uk.
 //
 //
 // A copy of the GNU General Public License is provided in the file gpl.txt. You
 // may also obtain a copy of the GNU General Public License on the World Wide
 // Web at http://www.gnu.org/licenses/gpl.html.
 */
package chemblast;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import jblast.chem.FormatFingerprintDatabase;
import jblast.chem.Mol2Sequence;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import jblast.sequence.BLAST;
import jblast.chem.FingerprintSequence;
import jblast.sequence.ProteinSequence;

/**
 *
 * @author Syed Asad Rahman, e-mail: s9asad@gmail.com
 */
public class ChemBlast {

    final static SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    final static String newline = System.getProperty("line.separator");

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        ChemBlast chemBlast = new ChemBlast();
        List<String> list = Arrays.asList(args);
        if (list.contains("-formatDB") && args.length == 3) {
            if (list.contains("-DB")) {
                int index = list.indexOf("-DB") + 1;
                String dbName = list.get(index);
                File f = new File(dbName);
                if (!f.exists()) {
                    System.err.println("ERROR: The DB name " + f.getAbsolutePath() + " not found!");
                } else {
                    System.out.println("INFO: The DB name is " + f.getAbsolutePath());
                    chemBlast.formatDB(f);
                }
            } else {
                System.err.println("WARNING: -DB option missing! ");
                printHelp();
                System.exit(1);
            }
        } else if (list.contains("-searchDB") && (args.length == 4 || args.length == 6)) {
            int index = list.indexOf("-searchDB") + 1;
            String querySMI = null;
            int topN = 10;
            if (list.contains("-querySMI")) {
                int smiIndex = list.indexOf("-querySMI") + 1;
                querySMI = list.get(smiIndex);
                if (list.contains("-top")) {
                    int topHitIndex = list.indexOf("-top") + 1;
                    topN = Integer.parseInt(list.get(topHitIndex).trim());
                }
            } else {
                printHelp();
                System.err.println("WARNING: -querySMI option missing! ");
                System.exit(1);
            }
            if (querySMI == null) {
                System.err.println("ERROR: Missing input SMILES!");
                printHelp();
                System.exit(1);
            }
            String dbName = list.get(index);
            File databasefile = new File(dbName);
            if (!databasefile.exists()) {
                System.err.println("ERROR: The DB name " + databasefile.getAbsolutePath() + " not found!");
            } else {
                System.out.println("INFO: The DB name is " + databasefile.getAbsolutePath());
                String fileName = databasefile.getName().split("\\.")[0];
                File indexfile = new File(databasefile.getParentFile(), fileName + ".idx");
                File formatDB = new File(databasefile.getParentFile(), fileName + ".fmt");
                if (!formatDB.exists() || !indexfile.exists()) {
                    chemBlast.formatDB(databasefile);
                } else {
                    System.out.println("INFO: Index file found " + indexfile.getAbsolutePath());
                    System.out.println("INFO: Formatted DB found " + formatDB.getAbsolutePath());
                }

                System.out.println("INFO: Reporting top " + topN + " hits for " + querySMI);
                /*
                 Perform BLAST search
                 */
                try {
//                    chemBlast.search(querySMI, indexfile, formatDB, topN);
                    chemBlast.searchParallel(querySMI, indexfile, formatDB, topN);
                } catch (InvalidSmilesException | IOException ex) {
                    Logger.getLogger(ChemBlast.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        } else {
            printHelp();
            System.exit(1);
        }
    }

    private void searchParallel(String querySmile, File indexfile, File formatDB, int top_hits)
            throws InvalidSmilesException, IOException {
        System.out.println("INFO: Searching DB...");
        IAtomContainer parseSmiles = smilesParser.parseSmiles(querySmile);
        FingerprintSequence query = Mol2Sequence.convertMol2Sequence("Query: " + querySmile, parseSmiles);
        ProteinSequence proteinSequence = new ProteinSequence(query.description(), query.elementsToString());
        try {
            System.out.println("INFO: Searching input database....");
            long t1 = System.currentTimeMillis();
            BLAST blast = new BLAST(proteinSequence, formatDB, indexfile, top_hits);
            long t2 = (System.currentTimeMillis());
            System.out.printf("Running Time: %.2f sec.", ((float) (t2 - t1) / (float) 1000));
            System.out.println("");
            System.out.println("INFO: Done");
        } catch (Exception ex) {
            Logger.getLogger(ChemBlast.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static void printHelp() {
        System.out.println("HELP!");
        System.out.println("1) Formatting the input database. Each line in the input file contains identifier for the molecule/drug and its smiles (Name\tSMILES).");
        System.out.println("\tjava -jar -formatDB -DB smiles.txt");
        System.out.println("2) Searching the input database. Input smiles file and the query smiles. This will return max. of 10 hits");
        System.out.println("\tjava -jar -searchDB smiles.txt -querySMI \"N1=C(C#N)N=C2C(=C1NCC(CC)C)NCN2CC3=CC=CC=C3\"");
    }

    public ChemBlast() {
    }
    /*
     Format the input DB- a) generate SMILES and b) transform it into tuples
     */

    void formatDB(File databasefile) {
        try {
            System.out.println("INFO: Formatting input database....");
            String fileName = databasefile.getName().split("\\.")[0];
            File indexfile = new File(databasefile.getParentFile(), fileName + ".idx");
            File formatDB = new File(databasefile.getParentFile(), fileName + ".fmt");
            String numberOfSeq = "-1";
            long t1 = System.currentTimeMillis();
            FormatFingerprintDatabase fingerprintDatabaseGenerator
                    = new FormatFingerprintDatabase(databasefile, indexfile, formatDB, numberOfSeq);
            long t2 = (System.currentTimeMillis());
            System.out.printf("Running Time: %.2f sec.", ((float) (t2 - t1) / (float) 1000));
            System.out.println("");
            System.out.println("\nINFO: Index file created " + indexfile.getAbsolutePath());
            System.out.println("INFO: DB formatted " + formatDB.getAbsolutePath());
            System.out.println("INFO: Done");
        } catch (IOException ex) {
            Logger.getLogger(ChemBlast.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
