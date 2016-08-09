#./BaMMmotif /home/wanwan/benchmark/JunD/order2/BaMM_new_fast_EM_changed /home/wanwan/benchmark/JunD/JunD.fasta --bindingSiteFile /home/wanwan/benchmark/JunD/motif_1.blocks --verbose --maxEMIterations 1 #--setSlow
./BaMMmotif /home/wanwan/benchmark/Bclaf1/order2_normal/BaMM_new_fast_EM_changed_alpha_vary /home/wanwan/benchmark/Bclaf1/Bclaf1.fasta --bindingSiteFile /home/wanwan/benchmark/Bclaf1/motif_1.blocks --verbose #--setSlow
#valgrind --leak-check=full --show-leak-kinds=all ./BaMMmotif /home/wanwan/benchmark/JunD/order2/BaMM_new_fast_r_changed /home/wanwan/benchmark/JunD/JunD.fasta --bindingSiteFile /home/wanwan/benchmark/JunD/motif_1.blocks --verbose

