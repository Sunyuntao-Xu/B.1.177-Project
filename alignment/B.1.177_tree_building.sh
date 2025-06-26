& "C:\iqtree-2.4.0-Windows\bin\iqtree2.exe" `
-s "E:\key_site_seqkit_filtered_B.1.177.countrydate_filtered.fasta" `
-m GTR+F+I+R3 `
-bb 1000 `
--bnni `
-nt AUTO `
-mem 0.95 `
-safe `
-pre B117_mainTree


/home/yuntao/Desktop/iqtree-3.0.1-Linux/bin/iqtree3_intel \
  -s "/media/yuntao/Toby Drive/key_site_seqkit_filtered_B.1.177.dedup.fasta" \
  -m GTR+F+I+R3 \
  --bnni \
  -nt AUTO \
  -mem 0.95 \
  -fast \
  -keep-ident \
  -pre /media/yuntao/Toby\ Drive/IQTREE_B.1.177/B117_mainTree

