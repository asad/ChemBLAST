# Usage
#sh format.sh data/database.txt "C[C@H](O)C(O)=O 5"
#

java -jar bin/ChemBlast.jar -searchDB $1 -querySMI $2 -top $3
