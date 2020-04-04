class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    counts = new int[codons.length];

    incCodons(inCodon);
    next = null;
  }

  private void incCodons(String s){
    for (int i = 0; i < codons.length; i++){
      if(codons[i].equals(s)){
        counts[i]++;
      }
    }

  }
  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops,
   * if not passes the task to the next node.
   * If there is no next node, add a new node to the list that would contain the codon.
   */
  private void addCodon(String inCodon){
    if(aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)){
      incCodons(inCodon);
    }

    else {
      if (next != null) {
        next.addCodon(inCodon);
      }
      else{
        next = new AminoAcidLL(inCodon);

      }
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for(int i = 0; i < counts.length; i++){
      sum+= counts[i];
    }
    return sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  public int sum(AminoAcidLL a){
    AminoAcidLL temp = a;
    int sum = 0;

    while(temp  != null){
      sum += temp.totalCount();
      temp = temp.next;
    }

    return sum;
  }
  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts.
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    if(inList.next == null && next == null){
      if(aminoAcid == inList.aminoAcid) {
        return totalDiff(inList);
      }
      return inList.totalCount() + totalCount();
    }

    if(aminoAcid == inList.aminoAcid){

      if(next == null) {
        return totalDiff(inList) + sum(inList.next);
      }

      if(inList.next == null) {
        return totalDiff(inList) + sum(next);
      }

      return totalDiff(inList) + next.aminoAcidCompare(inList.next);
    }

    if(this.next != null && aminoAcid < inList.aminoAcid){

      return totalCount() + next.aminoAcidCompare(inList);
    }

    return inList.totalCount() + aminoAcidCompare(inList.next);
  }

  /********************************************************************************************/
  /* Same as above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    if(inList.next == null && next == null){
      if(aminoAcid == inList.aminoAcid) {
        return totalDiff(inList);
      }
      return inList.totalCount() + totalCount();
    }
    if(aminoAcid == inList.aminoAcid){

      if(next == null) {
        return codonDiff(inList) + sum(inList.next);
      }

      if(inList.next == null) {
        return codonDiff(inList) + sum(next);
      }

      return codonDiff(inList) + next.codonCompare(inList.next);
    }
    if(next != null && aminoAcid < inList.aminoAcid){

      return totalCount() + next.codonCompare(inList);
    }

    return inList.totalCount() + codonCompare(inList.next);
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    // base case
    if(next == null){
      return new char[]{aminoAcid};
    }

    char[] a = next.aminoAcidList();
    char[] b = new char[a.length + 1];
    b[0] = aminoAcid;
    for (int i = 0; i < a.length; i++) {
      b[i+1] = a[i];
    }

    return b;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    // base case
    if(next == null){
      return new int[]{totalCount()};
    }

    int[] a = next.aminoAcidCounts();
    int[] ret = new int[a.length+1];
    ret[0] = totalCount();
    for(int i = 0; i < a.length; i++){
      ret[i+1] = a[i];
    }

    return ret;
  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    if(next == null){
      return true;
    }

    else if(aminoAcid > next.aminoAcid){
      return false;
    }

    return next.isSorted();
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    AminoAcidLL list = new AminoAcidLL(inSequence.substring(0,3));
    while(inSequence.length() > 3 && AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0,3)) != '*'){
      inSequence = inSequence.substring(3);
      list.addCodon(inSequence.substring(0,3));
    }
    return list;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    if(inList.next == null || inList == null){
      return inList;
    }

    char temp;
    for( AminoAcidLL i = inList; i.next != null; i = i.next){
      for (AminoAcidLL j = i.next; j != null; j = j.next){
        if(i.aminoAcid > j.aminoAcid){
          temp = i.aminoAcid;
          i.aminoAcid = j.aminoAcid;
          j.aminoAcid = temp;

        }
      }
    }
    return inList;
  }
}