#!/bin/bash

for cancer_type in $(cat /home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer_names.txt)
do
	cd /home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/$cancer_type/new/
	rsync -av --progress ./ feduati@biosim3.bmt.tue.nl:/home/feduati/Oscar/PhD_Oscar/new/data/PanCancer/$cancer_type/
done

