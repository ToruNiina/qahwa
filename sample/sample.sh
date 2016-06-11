#!/bin/sh

execute()
{
    echo "=======================================" >> ${logfilename}
    echo "execute: ./sample ${anchor_pos} ${coef}" >> ${logfilename}
    ./sample ${anchor_pos} ${coef} >> ${logfilename}
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> ${logfilename}
    echo "" >> ${logfilename}

    mv "sample_dist_${anchor_pos}_coef_${coef}.dcd" ./result/
}

main()
{
    : "if there is no ./result/ directory, make the directory" &&
    {
        if [ ! -d "result" ]; then
            mkdir result
        fi
    }

    : "prepair log file" && 
    {
        logfilename="sample.sh.log"
        echo "run sample" > ${logfilename}
    }

    : "run sample for some anch_pos" && 
    {
        for i in `seq 5.0 0.5 10.0`; do
            anchor_pos=${i}
            coef=10.0
            execute
        done
    }
    echo "end."
}

main $@
