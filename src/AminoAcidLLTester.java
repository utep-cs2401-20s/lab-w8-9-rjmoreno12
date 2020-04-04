import org.junit.jupiter.api.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;

public class AminoAcidLLTester {


    @Test
    // Tests aminoAcidCompare method
    public void aminoAcidCompareT1(){
        AminoAcidLL x = AminoAcidLL.createFromRNASequence("GAGGAGACCACCUGCGACUAG");
        AminoAcidLL y = AminoAcidLL.createFromRNASequence("GGUGGUGAGGAGGAGACCACCUAG");
        x = AminoAcidLL.sort(x);
        y = AminoAcidLL.sort(y);
        assertEquals(7, x.aminoAcidCompare(y));
    }

    @Test
    // Tests aminoAcidList method
    public void aminoAcidListT1(){
        char[] expected = {'F','L','P','G','K'};
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUUUUGCCCGGGAAA");
        assertArrayEquals(expected, test.aminoAcidList());
    }

    @Test
    // Tests aminoAcidList method COULD NOT GET IT TO WORK
    public void aminoAcidListT2(){
        char[] expected = {'V','M','I','K',};
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GUGAUGAUUAAG");
        assertArrayEquals(expected, test.aminoAcidList());
    }

    @Test
    // Tests aminoAcidCounts method
    public void aminoAcidCountsT1(){
        int[] expected = {3, 3, 3};
        String testSequence = "GCGGCGGCGUGUUGUUGUAAAAAAAAA";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence(testSequence);
        assertArrayEquals(expected, test.aminoAcidCounts());
    }

    @Test
    // Tests aminoAcidCounts method
    public void aminoAcidCountsT2(){
        int[] expected = {1, 1, 1};
        String testSequence = "GGAGAAGUC";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence(testSequence);
        assertArrayEquals(expected, test.aminoAcidCounts());
    }

    @Test
    // Tests sort method
    public void sortT1(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUCUUCUUAUUGCUUAUGAUACAGCAAUAA");
        test = AminoAcidLL.sort(test);
        assertEquals(true, test.isSorted());

    }

    @Test
    // Tests isSorted method with an unsorted list
    public void isSortedT1(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GUGGUAAUGAUCAUCCUCUUCUUCUCGUAA");
        assertEquals(false, test.isSorted());
    }

    @Test
    // Tests isSorted method with a sorted list
    public void isSortedT2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCGGCAUGUGACGAUGAAUUUGGGGGC");
        assertEquals(true, test.isSorted());
    }


    @Test
    // Tests createFromRNASequence with a stop codon in the middle of the sequence
    public void createFromRNASequenceT2(){
        String expected = "PQR";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("CCACAGCGUUAGGUGUGGAGCACGUAG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), test.aminoAcid);
            test = test.next;
        }
    }

    @Test
    public void sortT2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUC");
        test = AminoAcidLL.sort(test);
        assertEquals(true, test.isSorted());

    }

    @Test
    // Tests createFromRNASequence
    public void createFromRNASequenceT1(){
        String expected = "PQRVWST";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("CCCCAGCGUGUGUGGAGCACGUAG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), test.aminoAcid);
            test = test.next;
        }
    }

}