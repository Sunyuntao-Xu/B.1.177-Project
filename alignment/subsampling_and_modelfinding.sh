seqkit sample -n 5000 -s 42 E:/key_site_seqkit_filtered_B.1.177.dedup.fasta -o E:/subsamples/B117_dedup_5k.fasta

& "C:\iqtree-2.4.0-Windows\bin\iqtree2.exe" `
  -s E:\subsamples\B117_dedup_5k.fasta `
  -m MFP `
  -nt AUTO `
  -pre E:\modeltest\B117_5k_modeltest
