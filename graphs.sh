#python3 plot.py output.txt "Simulated Medulloblastoma Sample" "graphs/sms.jpeg"
python3 plot.py output.txt "Simulated Medulloblastoma Sample (100bp mean)" "graphs/sms100mean.jpeg" -n 100 mean
python3 plot.py output.txt "Simulated Medulloblastoma Sample (100bp median)" "graphs/sms100med.jpeg" -n 100 median
python3 plot.py output.txt "Simulated Medulloblastoma Sample (100bp 95% quantile)" "graphs/sms100q95.jpeg" -n 100 quantile .95
python3 plot.py output.txt "Simulated Medulloblastoma Sample (10kbp mean)" "graphs/sms10kmean.jpeg" -n 10000 mean

#python3 plot.py output.txt "Control Sample" "graphs/cs.jpeg"
python3 plot.py results-wt/output.txt "Control Sample (100bp mean)" "graphs/cs100mean.jpeg" -n 100 mean
python3 plot.py results-wt/output.txt "Control Sample (100bp median)" "graphs/cs100med.jpeg" -n 100 median
python3 plot.py results-wt/output.txt "Control Sample (100bp 95% quantile)" "graphs/cs100q95.jpeg" -n 100 quantile .95
python3 plot.py results-wt/output.txt "Control Sample (10kbp mean)" "graphs/cs10kmean.jpeg" -n 10000 mean
