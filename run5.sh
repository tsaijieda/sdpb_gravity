#list=(0.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0)
#list=(1.7)
#list=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 3.0 5.0 7.0 9.0 11.0 13.0)
list=(1.0 1.2 1.5 1.7 2.0 2.2 2.5 2.7 3.0 3.5 4.0)
#list=(1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0)
#list=(0.0 10.0 50.0 100.0 -10.0 -50.0 -100.0 200.0 -200.0)
#list=(2000.0 1000.0 500.0)
#list=(1.0 4.5 -4.0 6.0 8.0 10.0 -6.0 -8.0 -10.0 20.0 50.0 -20.0 -50.0)
#list=(2.0 -7.0)
#this test starts from here
#list=(-200.0 -100.0 -50.0 -20.0 0.0 20.0 50.0 100.0 200.0 400.0)
#list=(400.0 1000.0 1500.0 2000.0)
#list=(10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40)
#list=(-1.0 -2.0 -3.0 -4.0 -5.02)
#list=(0.0 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0)
#list=(-1.5 -1.4 -1.3 -1.2 -1.1)
#list=(0.25 0.50 0.75 1.0 1.25 1.50 1.75 2.0 2.25 2.50 2.75 3.0 3.25 3.50 3.75 4.0 4.25 4.50 4.75 5.0 5.25 5.50 5.75 6.0 6.25 6.50 6.75 7.0)
#list=(-0.75 -0.25)
spinlist=(0.0)
dir=$(pwd)
echo "$dir"
#rm -rf inter_result.txt
for ((j = 0; j < ${#spinlist[@]}; j++))
do
for ((i = 0; i < ${#list[@]}; i++)) 
    do
        echo "Start process"
	rm -r out.ck
        #rm -rf  sdp.zip.ck  sdp.zip_out  sdp.zip  division.txt result.txt

        mu=${list[i]}
	spin0=${spinlist[j]}

        python3 change.py range.txt $mu -1 $spin0

        echo "Generating xml file"
	wolframscript -file "gravity-supergrav-improved-cont.m"

        echo "Turning xml to sdp file"
        docker run --shm-size=8192m -v "$dir":/data bootstrapcollaboration/sdpb:master pmp2sdp 2048 -i /data/out.json -o /data/out

	echo "Running sdpb"
        docker run --shm-size=8192m -v "$dir":/data bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 48 sdpb --writeSolution="x,y,z,X,Y" --maxComplementarity=1e1000 --dualityGapThreshold=1e-30 --primalErrorThreshold=1e-30 --dualErrorThreshold=1e-30 --precision=2048 --procsPerNode=32 --maxIterations=50000 -s /data/out
        
	python3 store_result.py out_out/out.txt "result/sarahtest.txt"

        #rm -rf upper_spectrum 

        #cp -rf sdp.zip_out upper_spectrum  

        #rm -rf sdp.zip.ck  sdp.zip_out  sdp.zip
    done
done

