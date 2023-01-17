function accept(rec) {
  //Only mapped reads were kept
  if (rec.getReadUnmappedFlag()) return false;
  //Mapping quality should be larger than 0
  //MPXV specific ITR (bwa mapping quality equal 0)
  if (rec.getMappingQuality() < 0) return false;
  //aligned proportion > 80%
  numAlignedBases = 0;
  b = rec.getAlignmentBlocks();
  //print(b.size())
  for (i=0 ; i < b.size(); i++) {
     //a =  alignmentBlock.getLength();
     //numAlignedBases += a;
     numAlignedBases += b[i].getLength();
     //print( rec.getReadName()  )
  }
  //print(numAlignedBases / rec.getReadLength());
  return  numAlignedBases / rec.getReadLength() > 0.8;
}
accept(record);
