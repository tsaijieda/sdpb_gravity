	dir=$(pwd)
	echo "Start process"
	rm -r out.ck
        #rm -rf  sdp.zip.ck  sdp.zip_out  sdp.zip  division.txt result.txt

        echo "Generating xml file"
        wolframscript -file gravity.m

        echo "Turning xml to sdp file"
        docker run --shm-size=4096m -v "$dir":/data bootstrapcollaboration/sdpb:master pmp2sdp 1024 -i /data/out.json -o /data/out

	echo "Running sdpb"
        docker run --shm-size=4096m -v "$dir":/data bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 48 sdpb --precision=1024 --procsPerNode=32 --maxIterations=5000 -s /data/out

        #rm -rf upper_spectrum 

        #cp -rf sdp.zip_out upper_spectrum  

        #rm -rf sdp.zip.ck  sdp.zip_out  sdp.zip

