nohup docker run --rm --gpus all \
-v $HOME/models:$HOME/models \
-v $HOME/outputs:$HOME/outputs \
-v /batch-output/working_dir:/input \
rfdiffusion \
inference.output_prefix=/input//RFDiffusion_prompt1 \
inference.model_directory_path=$HOME/models \
inference.ckpt_override_path=$HOME/models/ActiveSite_ckpt.pt \
inference.input_pdb=/input/RFDiffusion/antibody_antigen_crystalStructure.pdb \
inference.num_designs=200 \
'contigmap.contigs=[1-100/A235-236/1-100/A307-309/1-100/A478-483/1-100/A510-511/1-100]' \
'contigmap.length=60-700'  >/batch-output/working_dir/RFDiffusion/prompt1.out &

nohup docker run --rm --gpus all \
-v $HOME/models:$HOME/models \
-v $HOME/outputs:$HOME/outputs \
-v /batch-output/working_dir:/input \
rfdiffusion \
inference.output_prefix=/input//RFDiffusion_prompt2 \
inference.model_directory_path=$HOME/models \
inference.ckpt_override_path=$HOME/models/ActiveSite_ckpt.pt \
inference.input_pdb=/input/RFDiffusion/antibody_antigen_crystalStructure.pdb \
inference.num_designs=200 \
'contigmap.contigs=[1-100/A235-236/1-100/A307-309/1-100/A478-484/1-100/A509-513/1-100]' \
'contigmap.length=60-700'  >/batch-output/working_dir/RFDiffusion/prompt2.out &

##calculate radius of gyration
python3 /batch-output/202409_RFDiffusion/scripts/radiusOfGyration_folder.py antigen/RFDiffusion/prompt2
python3 /batch-output/202409_RFDiffusion/scripts/radiusOfGyration_folder.py antigen/RFDiffusion/prompt2
###change the jsonl file to do proteinMPNN on the fixed positions
python3 /batch-output/202409_RFDiffusion/scripts/processTRBinAFolder.py ./antigen/RFDiffusion/prompt1 ./antigen/RFDiffusion/prompt1TRBout.txt
###measure the distance between the proteinMPNN generated structures and the original, use that as filter??

python ./scripts/calculateDistanceBetweenMoftifs.py ./antigen/RFDiffusion/prompt2_proteinMotifPositions.txt ./antigen/RFDiffusion/prompt2MotifDistance.txt
##summarize lengths of the proteins generated

python ./scripts/calculateDistanceBetweenMoftifs.py ./antigen/RFDiffusion/prompt2_proteinMotifPositions.txt 

sed -i -e 's/":/.pdb":/g' ./antigen/RFDiffusion/prompt2TRBout.txt
python ./scripts/calculateDistanceBetweenMoftifs.py ./antigen/RFDiffusion/prompt2TRBout.txt ./antigen/RFDiffusion/prompt2Distances.txt


##alphafold the ones that we do like (ie good radius of gyration and very little distance to the original, based only on RFDiffusion)
##use alphafold to check the confidence in the proteinMPNN predictions?
#(and use plddt as a filter for the ones that we're confident in?)

###summarize difference between the different proteinMPNN structures by calcualting the distance between 2 PDB structures


##run 100 on the first prompt? see the range of scores we get and also make sure we can do it on the variable lentsh

path_for_parsed_chains="/output/ProteinMPNN/parsed_pdbsPrompt1All.jsonl"
path_for_assigned_chains="/output/ProteinMPNN/assigned_pdbsPrompt1All.jsonl"
path_for_fixed_positions="/output/ProteinMPNN/fixed_pdbsPrompt1All.jsonl"
folder_with_pdbs="/output/RFDiffusion/prompt1HigherThan1"
chains_to_design="A"

docker run --gpus all \
-v /batch-output/working_dir//:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
--input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
--input_path=$path_for_parsed_chains \
--output_path=$path_for_assigned_chains --chain_list "$chains_to_design"


nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/ProteinMPNN/prompt1All.jsonl" \
--out_folder "/output/ProteinMPNN/prompt1AllHigher" \
--num_seq_per_target 100 \
--sampling_temp "0.1" \
--seed 37 \
--batch_size 1 >/batch-output/working_dir//ProteinMPNN/proteinMPNNPrompt1AllHigher.txt &


path_for_parsed_chains="/output/ProteinMPNN/parsed_pdbsPrompt2All.jsonl"
path_for_assigned_chains="/output/ProteinMPNN/assigned_pdbsPrompt2All.jsonl"
#path_for_fixed_positions="/output/ProteinMPNN/fixed_pdbsPrompt2All.jsonl"
folder_with_pdbs="/output/RFDiffusion/prompt2"
chains_to_design="A"

docker run --gpus all \
-v /batch-output/working_dir//:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
--input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
--input_path=$path_for_parsed_chains \
--output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/ProteinMPNN/prompt2TRBout.jsonl" \
--out_folder "/output/ProteinMPNN/prompt2" \
--num_seq_per_target 100 \
--sampling_temp "0.1" \
--seed 37 \
--batch_size 1 >/batch-output/working_dir/ProteinMPNN/proteinMPNNPrompt2All.txt &


###
 nohup colabfold_batch --model-type alphafold2 \
 /batch-output/working_dir/ProteinMPNN/prompt1All/seqs/RFDiffusion_prompt1_1.fa \
 /batch-output/working_dir/ProteinMPNN/Prompt1AllFolding/RFDiffusion_prompt1_1_folding \
  >/batch-output/working_dir//ProteinMPNN/Prompt1AllFolding.txt &



 nohup colabfold_batch --model-type alphafold2 \
 /batch-output/working_dir/ProteinMPNN/prompt1All/seqs/RFDiffusion_prompt1_199.fa \
 /batch-output/working_dir/ProteinMPNN/Prompt1AllFolding/RFDiffusion_prompt1_199_folding \
  >/batch-output/working_dir//ProteinMPNN/Prompt1AllFolding199.txt &
  ##very bad plddt for prompt1_199
  ##we shoudl check the scores? as well check whether higher temperature helps?
  
  
  nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/ProteinMPNN/prompt1All.jsonl" \
--out_folder "/output/ProteinMPNN/prompt1All2" \
--num_seq_per_target 100 \
--sampling_temp "0.2" \
--seed 37 \
--batch_size 10 >/batch-output/working_dir//ProteinMPNN/proteinMPNNPrompt1All2.txt &

##try prompt1_199 temperature0.2?

nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/ProteinMPNN/prompt1All.jsonl" \
--out_folder "/output/ProteinMPNN/prompt1All2" \
--num_seq_per_target 100 \
--sampling_temp "0.2" \
--seed 37 \
--batch_size 10 >/batch-output/working_dir//ProteinMPNN/proteinMPNNPrompt1All2.txt &


python ./scripts/outputDistanceFromPDBandMofit2.py ./antigen/RFDiffusion/prompt2_proteinMotifPositions.csv ./antigen/RFDiffusion/prompt2_proteinMotifPositions_output2.csv
 
 
##use the above script on the original pdb
python ./scripts/outputDistanceFromPDBandMofit2.py ./antigen/RFDiffusion/'prompt2_proteinMotifPositions=3.csv' ./antigen/originalantigenAndMotifPositionOutput.csv ##i think this is prompt1
 
python ./scripts/outputDistanceFromPDBandMofit2.py ./antigen//originalantigenAndMotifPosition2.csv ./antigen/originalantigenAndMotifPositionOutput2.csv
##find the difference between these 2 matrices
#(I guess i could actually use R for it but we'll just continue chatgpt-ing)
python ./scripts/distanceBetween2pdb.py ./antigen/RFDiffusion/prompt2_proteinMotifPositions_output2.csv ./antigen/originalantigenAndMotifPositionOutput2.csv ./antigen/output_differences2.csv
##add 2 for the ord2 (ie supposedly eucliean distance)


##doing the same for prompt1
python ./scripts/outputDistanceFromPDBandMofit2.py ./antigen/RFDiffusion/prompt1_motifPositions ./antigen/RFDiffusion/prompt1_proteinMotifPositions_output.csv
#python ./scripts/outputDistanceFromPDBandMofit2.py ./antigen//originalantigenAndMotifPosition2.csv ./antigen/originalantigenAndMotifPositionOutput2.csv
python ./scripts/distanceBetween2pdb.py ./antigen/RFDiffusion/prompt1_proteinMotifPositions_output.csv ./antigen/originalantigenAndMotifPositionOutput.csv ./antigen/q.csv


##prompt3 with NERK
nohup docker run --rm --gpus all \
-v $HOME/models:$HOME/models \
-v $HOME/outputs:$HOME/outputs \
-v /batch-output/working_dir:/input \
rfdiffusion \
inference.output_prefix=/input//RFDiffusion_prompt3 \
inference.model_directory_path=$HOME/models \
inference.ckpt_override_path=$HOME/models/ActiveSite_ckpt.pt \
inference.input_pdb=/input/RFDiffusion/antibody_antigen_crystalStructure.pdb \
inference.num_designs=200 \
'contigmap.contigs=[1-100/A234-236/1-100/A307-310/1-100/A315-317/1-100/A478-484/1-100/A492-494/1-100/A510-511/1-100]' \
'contigmap.length=60-700'  >/batch-output/working_dir/RFDiffusion/prompt3.out &

nohup docker run --rm --gpus all \
-v $HOME/models:$HOME/models \
-v $HOME/outputs:$HOME/outputs \
-v /batch-output/working_dir:/input \
rfdiffusion \
inference.output_prefix=/input//RFDiffusion_prompt2 \
inference.model_directory_path=$HOME/models \
inference.ckpt_override_path=$HOME/models/ActiveSite_ckpt.pt \
inference.input_pdb=/input/RFDiffusion/antibody_antigen_crystalStructure.pdb \
inference.num_designs=200 \
'contigmap.contigs=[1-100/A235-236/1-100/A307-309/1-100/A478-484/1-100/A509-513/1-100]' \
'contigmap.length=60-700'  >/batch-output/working_dir/RFDiffusion/prompt2.out &




####testing the limit of proteinmppnn
##prompt1, trying cutoff 1.09ish
#RFDiffusion_prompt1_169

 nohup colabfold_batch --model-type alphafold2 \
 /batch-output/working_dir/ProteinMPNN/prompt1All/seqs/RFDiffusion_prompt1_169.fa \
 /batch-output/working_dir/ProteinMPNN/Prompt1AllFolding/RFDiffusion_prompt1_169_folding \
  >/batch-output/working_dir//ProteinMPNN/Prompt1AllFolding169.txt &
  
##may consider only using the best scores (ie lowest) for all 200
##i think we should have a list of the pdb files generated from RFDiffusion, and compare it to the best scores structures
##from proteinMPNN, and see if they're similar in structure, at least at the Calpha
##and maybe convert from cif to pdb format 


##running proteinMPNN on prompt3
path_for_parsed_chains="/output/ProteinMPNN/parsed_pdbsPrompt3All.jsonl"
path_for_assigned_chains="/output/ProteinMPNN/assigned_pdbsPrompt3All.jsonl"
folder_with_pdbs="/output/RFDiffusion/prompt3"
chains_to_design="A"

docker run --gpus all \
-v /batch-output/working_dir//:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
--input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
--input_path=$path_for_parsed_chains \
--output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python3 ./scripts/processTRBinAFolder.py ./antigen/RFDiffusion/prompt3 ./antigen/RFDiffusion/prompt3TRBout.txt


nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/RFDiffusion/prompt3TRBout.jsonl" \
--out_folder "/output/ProteinMPNN/prompt3" \
--num_seq_per_target 100 \
--sampling_temp "0.1" \
--seed 37 \
--batch_size 1 >/batch-output/working_dir/ProteinMPNN/proteinMPNNPrompt3All.txt &


##checking the score for all prompt1
##alphafold
nohup colabfold_batch --model-type alphafold2 \
/batch-output/working_dir/AF/prompt1FastaToRunAFon.fa \
/batch-output/working_dir/AF/prompt1FastaToRunAFon \
>/batch-output/working_dir/AF/prompt1FastaToRunAFon.txt &
  
  
##AF3 on the G6E server
docker run -d -v /batch-output/202409_RFDiffusion/AF3selfRun:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--json_path=/af_input/fold_prompt1Test.json.json --model_dir=/root/models \
--output_dir=/af_input/fold_prompt1Test

###checking against iedb
##ie grabbing the sequences to blast again
#go through the xml files in /batch-output/database
nohup colabfold_batch --model-type alphafold2 \
/batch-output/working_dir/AF/prompt2FastaToRunAFon.fa \
/batch-output/working_dir/AF/prompt2/prompt2FastaToRunAFon \
>/batch-output/working_dir/AF/prompt2FastaToRunAFon.txt &

##prompt3
nohup colabfold_batch --model-type alphafold2 \
/batch-output/working_dir/AF/prompt3FastaToRunAFon.fa \
/batch-output/working_dir/AF/prompt3/prompt3FastaToRunAFon \
>/batch-output/working_dir/AF/prompt3FastaToRunAFon.txt &
  


docker run -v /batch-output/database/:/database/ \
620901718958.dkr.ecr.us-east-2.amazonaws.com/blast makeblastdb \
-dbtype 'prot' -in /database/GenbankIDAndSequence2.fa 



  
###

for file in ls *.xml 
do
echo $file >>../GenbankIDAndSequenceAll2.fa
grep -e '<GenBankId>' -e '<LinearSequence>' $file >>../GenbankIDAndSequenceAll2.fa
done

python /batch-output/202409_RFDiffusion/scripts/parseXML4.py \
/batch-output/database/GenbankIDAndSequenceAll2.fa \
/batch-output/database/GenbankIDAndSequenceAllDB.fa

nohup docker run -v /batch-output/database/:/database/ \
620901718958.dkr.ecr.us-east-2.amazonaws.com/blast makeblastdb \
-dbtype 'prot' -in /database/GenbankIDAndSequenceAllDB.fa &
###blastp against the current database

docker run -v /batch-output/database/:/database/ \
-v /batch-output/working_dir/AF/:/input/ \
620901718958.dkr.ecr.us-east-2.amazonaws.com/blast blastp \
-db /database/GenbankIDAndSequenceAllDB.fa -outfmt 6 -num_threads 6 \
-query /input/prompt2FastaToRunAFon.fa -evalue 0.01 \
>/batch-output/working_dir/AF/prompt2AgainstIEDB.txt



docker run -v /batch-output/database/:/database/ \
-v /batch-output/working_dir/AF/:/input/ \
620901718958.dkr.ecr.us-east-2.amazonaws.com/blast blastp \
-db /database/GenbankIDAndSequenceAllDB.fa -outfmt 6 -num_threads 6 \
-query /input/prompt2FastaToRunAFon.fa -evalue 0.1 \
>/batch-output/working_dir/AF/prompt2AgainstIEDB.txt

##checking that the AF2 output structure is the structure we think we want?
##i guess for prompt1 since 2 and 3 are not done yet

docker run -d -v /batch-output/202409_RFDiffusion/AF3selfRun:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--json_path=/af_input/fold_prompt1Test.json --model_dir=/root/models \
--output_dir=/af_input/fold_prompt1Test
 
 nohup colabfold_batch --model-type alphafold2 \
 /batch-output/working_dir/AF/prompt1Higher/prompt1HigherFastaToRunAFon.fa \
 /batch-output/working_dir/AF/prompt1Higher \
  >/batch-output/working_dir/AF/prompt1HigherFastaToRunAFon.txt &
  
##let's try to blast against iedb

docker run -v /batch-output/database/:/database/ \
-v /batch-output/working_dir/AF/:/input/ \
620901718958.dkr.ecr.us-east-2.amazonaws.com/blast blastp \
-db /database/GenbankIDAndSequenceAllDB.fa -outfmt 6 -num_threads 6 \
-query /input/prompt1Higher/prompt1HigherFastaToRunAFon.fa -evalue 0.1 \
>/batch-output/working_dir/AF/prompt1HigherAgainstIEDB.txt

nohup docker run -v /batch-output/working_dir/AF3:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--input_dir =/af_input/bestPLDDT10 --model_dir =/root/models \
--output_dir =/af_input/bestPLDDT10Output >bestPLDDT10.log &


nohup docker run -v /batch-output/working_dir/AF3:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--input_dir=/af_input/prompt2Jsons --model_dir=/root/models \
--output_dir=/af_input/bestPLDDTprompt2 >bestPLDDTprompt2.log &


for file in ls *json
  do
  echo $file
  sed -i -e 's/\[\n{/\[{/g' $file
done
  
  
  
##get the other fastas of proteins that look god
nohup docker run -v /batch-output/working_dir/AF3:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--input_dir=/af_input/prompt1GoodPAEJsons --model_dir=/root/models \
--output_dir=/af_input/bestPAEprompt1 >prompt1GoodPAEJsons.log &

nohup docker run -v /batch-output/working_dir/AF3:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--input_dir=/af_input/prompt1GoodPAEJsonsII --model_dir=/root/models \
--output_dir=/af_input/bestPAEprompt1II >prompt1GoodPAEJsonsII.log &


nohup docker run --rm --gpus all \
-v $HOME/models:$HOME/models \
-v $HOME/outputs:$HOME/outputs \
-v /batch-output/working_dir:/input \
rfdiffusion \
inference.output_prefix=/input//RFDiffusion_prompt1ii \
inference.model_directory_path=$HOME/models \
inference.ckpt_override_path=$HOME/models/ActiveSite_ckpt.pt \
inference.input_pdb=/input/RFDiffusion/antibody_antigen_crystalStructure.pdb \
inference.num_designs=400 \
'contigmap.contigs=[1-100/A235-236/1-100/A307-309/1-100/A478-483/1-100/A510-511/1-100]' \
'contigmap.length=60-700'  >/batch-output/working_dir/RFDiffusion/prompt1ii.out &

nohup docker run --rm --gpus all \
-v $HOME/models:$HOME/models \
-v $HOME/outputs:$HOME/outputs \
-v /batch-output/working_dir:/input \
rfdiffusion \
inference.output_prefix=/input//RFDiffusion_prompt2ii \
inference.model_directory_path=$HOME/models \
inference.ckpt_override_path=$HOME/models/ActiveSite_ckpt.pt \
inference.input_pdb=/input/RFDiffusion/antibody_antigen_crystalStructure.pdb \
inference.num_designs=400 \
'contigmap.contigs=[1-100/A235-236/1-100/A307-309/1-100/A478-484/1-100/A509-513/1-100]' \
'contigmap.length=60-700'  >/batch-output/working_dir/RFDiffusion/prompt2ii.out &

python3 /batch-output/202409_RFDiffusion/scripts/processTRBinAFolder.py \
./antigen/RFDiffusion/RFDiffusion_prompt1ii ./antigen/RFDiffusion/prompt1iiTRBout.txt


python3 /batch-output/202409_RFDiffusion/scripts/processTRBinAFolder.py \
./antigen/RFDiffusion/RFDiffusion_prompt2ii ./antigen/RFDiffusion/prompt2iiTRBout.txt

path_for_parsed_chains="/output/ProteinMPNN/parsed_pdbsPrompt1iiAll.jsonl"
path_for_assigned_chains="/output/ProteinMPNN/assigned_pdbsPrompt1iiAll.jsonl"
folder_with_pdbs="/output/RFDiffusion/RFDiffusion_prompt1ii"
chains_to_design="A"

docker run --gpus all \
-v /batch-output/working_dir//:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
--input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
--input_path=$path_for_parsed_chains \
--output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/ProteinMPNN/prompt1ii.jsonl" \
--out_folder "/output/ProteinMPNN/prompt1ii" \
--num_seq_per_target 100 \
--sampling_temp "0.1" \
--seed 37 \
--batch_size 1 >/batch-output/working_dir//ProteinMPNN/proteinMPNNPrompt1ii.txt &

path_for_parsed_chains="/output/ProteinMPNN/parsed_pdbsPrompt2iiAll.jsonl"
path_for_assigned_chains="/output/ProteinMPNN/assigned_pdbsPrompt2iiAll.jsonl"
folder_with_pdbs="/output/RFDiffusion/RFDiffusion_prompt2ii"
chains_to_design="A"

docker run --gpus all \
-v /batch-output/working_dir//:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
--input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
--input_path=$path_for_parsed_chains \
--output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

nohup docker run --gpus all \
-v /batch-output/working_dir/:/output \
protein_mpnn python /app/ProteinMPNN/protein_mpnn_run.py \
--jsonl_path $path_for_parsed_chains \
--chain_id_jsonl $path_for_assigned_chains \
--fixed_positions_jsonl "/output/ProteinMPNN/prompt2ii.jsonl" \
--out_folder "/output/ProteinMPNN/prompt2ii" \
--num_seq_per_target 100 \
--sampling_temp "0.1" \
--seed 37 \
--batch_size 1 >/batch-output/working_dir//ProteinMPNN/proteinMPNNPrompt2ii.txt &


##let's do solubility and then binding
#PLM-sol for solubility
conda activate /home/ubuntu/.conda/envs/PLM_Sol
cd embedding_datset
#Change the file path (sequences_file: ./Train_dataset.fasta prefix: ./Train_dataset_emb)
bio_embeddings embedding_protT5.yml
 
#haddock3 for binding


###prompt1ii
nohup docker run -v /batch-output/working_dir/AF3:/af_input \
-v /batch-output/202409_RFDiffusion/AF3selfRun/models:/root/models \
-v /home/ec2-user/databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py \
--input_dir=/af_input/prompt1iiProteinMPNNgoodJsons --model_dir=/root/models \
--output_dir=/af_input/bestPAEprompt1ii >prompt1iiProteinMPNNgoodJsons.log &

####
nohup colabfold_batch --model-type alphafold2 \
/batch-output/working_dir/ProteinMPNN/prompt1iiToRunAFon.fa \
/batch-output/working_dir/AF/prompt1iiFastaToRunAFon \
>/batch-output/working_dir/AF/prompt1iiFastaToRunAFon.txt &

nohup colabfold_batch --model-type alphafold2 \
/batch-output/working_dir/AF/prompt2iiToRunAFon.fa \
/batch-output/working_dir/AF/prompt2iiFastaToRunAFon \
>/batch-output/working_dir/AF/prompt2iiFastaToRunAFon.txt &



###running the haddock

prompt1iiToRunAFon.fa 

###running blastp
docker run -v /batch-output/database/:/database/ \
-v /batch-output/working_dir/best6/:/input/ \
620901718958.dkr.ecr.us-east-2.amazonaws.com/blast blastp \
-db /database/GenbankIDAndSequenceAllDB.fa -outfmt 6 -num_threads 6 \
-query /input/sixToTest.fasta -evalue 0.01 \
>/batch-output/working_dir/best6/best6AgainstIEDB.txt
