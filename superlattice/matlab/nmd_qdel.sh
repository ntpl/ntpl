MAXJOBID=578205
i="578004"
while [ $i -le $MAXJOBID ]
do
	qdel $i

i=$(( $i + 1 ))
done


