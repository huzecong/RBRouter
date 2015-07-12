rm testcase/*.ps
for FILE in `ls testcase`; do
	NAME=`echo $FILE | cut -d. -f1`;
	./RBRouter testcase/$FILE testcase/$NAME.ps;
done
