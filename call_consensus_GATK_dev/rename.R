require(Biostrings)
args = commandArgs(T)
fa = args[1]
rename = args[2]

fa_seq = readDNAStringSet(fa)
names(fa_seq) = rename

writeXStringSet(fa_seq, fa)
