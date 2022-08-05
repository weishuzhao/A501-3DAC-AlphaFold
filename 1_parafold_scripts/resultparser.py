import os

dirname = "3DAC"

os.chdir("3DAC_result")
os.mkdir("3DAC_summary")
flog = open("parser.log","w")

for i in range(1,303):
    os.mkdir("3DAC_summary/3DAC_%.3d"%i)
    os.system("cp 3DAC_%.3d/*.err 3DAC_summary/3DAC_%.3d" %(i,i))
    os.system("cp 3DAC_%.3d/*.out 3DAC_summary/3DAC_%.3d" %(i,i))
    if os.path.exists("3DAC_%.3d/ranking_debug.json"%i) == True:
        flog.write("3DAC_%.3d Finished\n"%i)
        os.system("cp 3DAC_%.3d/ranking_debug.json 3DAC_summary/3DAC_%.3d" %(i,i))
        os.system("cp 3DAC_%.3d/relaxed_model_1.pdb 3DAC_summary/3DAC_%.3d" %(i,i))
        os.system("cp 3DAC_%.3d/relaxed_model_2.pdb 3DAC_summary/3DAC_%.3d" %(i,i))
        # os.system("cp 3DAC_%.3d/result_model_1.pkl 3DAC_summary/3DAC_%.3d" %(i,i))
        # os.system("cp 3DAC_%.3d/result_model_2.pkl 3DAC_summary/3DAC_%.3d" %(i,i))
    else:
        flog.write("3DAC_%.3d Error\n"%i)


