from pymol import cmd
import os
import time

start_time = time.time()

index_total_1 = 1878
# index_list_1 = range(1,index_total_1+1)
index_list_1 = [1,60]
index_total_2 = 1196

rmsd_list_total = []
for index_1 in index_list_1:
    print("Start RMSD calculation of index 1: %.3d"%index_1)
    if os.path.exists("A501_non_result/A501_non_%.3d/unrelaxed_model_1.pdb"%index_1):
        f_out = open("RMSD/RMSD_A501_non_%.3d.txt"%index_1,"w")
        for index_2 in range(1,index_total_2+1):
            try:
                cmd.load("A501_non_result/A501_non_%.3d/unrelaxed_model_1.pdb"%index_1,"A501_1")
                cmd.load("3DAC_non_result/3DAC_non_%.3d/unrelaxed_model_1.pdb"%index_2,"3DAC_1")
                f_out.write("A501_%.3d,3DAC_%.3d"%(index_1,index_2)+str(cmd.align("A501_1","3DAC_1",cycles=5))+"\n")
                cmd.reinitialize()
                print("%.3d"%index_2,end=",")
            except:
                f_out.write("A501_%.3d,3DAC_%.3d"%(index_1,index_2)+"(No Data)"+"\n")
                print("Failed %.3d"%index_2,end=",")
        f_out.close()
        print("\nFinished\n")
    else:
        f_out = open("RMSD/RMSD_A501_non_%.3d.txt"%index_1,"w")
        f_out.close()
        print("\nFailed\n")

end_time = time.time() - start_time
print(end_time)