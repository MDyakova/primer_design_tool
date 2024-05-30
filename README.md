# primer_design_tool

docker build -t dockerfile:primerdesign ./ 
docker run -p 5000:5000 -d dockerfile:primerdesign

## docker hub
docker tag dockerfile:primerdesign mdyakova/primer_design_tool:v1
docker login
docker push mdyakova/primer_design_tool:v1

## Check code quality
python3 -m black main.py
python3 lint.py --threshold 7
