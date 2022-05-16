#准备FcC_supermatrix.fas文件，注意：里面的物种序列排列必须是按首字母顺序排，不然后面一旦删除种名将混乱
#将按首字母排的序列文件，删除种名（例如：删除 >Apis_xxxxx）,仅保留48条序列（这一步若物种不多、文件不大可以手动删除），保存为sequence_withoutname
#将sequence_withoutname文件每个字符中间加入一个空格，生成sequences.fina文件
cat sequence_withoutname | sed 's/.\{1\}/& /g' > sequences.fina

#将sequences.fina文件进行行转列，生成seq_shu文件
cat sequences.fina | awk '{$1="";print}' | awk 'BEGIN{c=0;} {for(i=1;i<=NF;i++) {num[c,i] = $i;} c++;} END{ for(i=1;i<=NF;i++){str=""; for(j=0;j<NR;j++){ if(j>0){str = str" "} str= str""num[j,i]}printf("%s\n", str)} }' > seq_shu

#将seq_shu文件随机取10000行（即随机取10000个位点），生成seq_cv文件
shuf -n10000 seq_shu > seq_cv
rm seq_shu

#将seq_cv文件中的10000行均分为10份
cat seq_cv | sed -n '1,1000p' > 01
cat seq_cv | sed -n '1001,2000p' > 02
cat seq_cv | sed -n '2001,3000p' > 03
cat seq_cv | sed -n '3001,4000p' > 04
cat seq_cv | sed -n '4001,5000p' > 05
cat seq_cv | sed -n '5001,6000p' > 06
cat seq_cv | sed -n '6001,7000p' > 07
cat seq_cv | sed -n '7001,8000p' > 08
cat seq_cv | sed -n '8001,9000p' > 09
cat seq_cv | sed -n '9001,10000p' > 10

#将01-10这10个分区文件，列转行，并删除空格
sed 's/[ ][ ]*/,/g' 01 > 0101
awk -F, '{for(i=1;i<=NF;i=i+1){a[NR,i]=$i}}END{for(j=1;j<=NF;j++){str=a[1,j];for(i=2;i<=NR;i++){str=str " " a[i,j]}print str}}' 0101 > 010101
sed 's/[ ][ ]*//g' 010101 > 01010101
rm 0101 010101

#将每行序列加上相对应的物种名（01010101手动替换）
sed -i '1i\>Acerentomon_sp' 01010101
sed -i '3i\>Acontista_multicolor' 01010101
sed -i '5i\>Anopheles_arabiensis' 01010101
sed -i '7i\>Apachyus_charteceus' 01010101
sed -i '9i\>Apis_mellifera' 01010101
sed -i '11i\>Atelura_formicaria' 01010101
sed -i '13i\>Bombyx_mori' 01010101
sed -i '15i\>Boreus_hyemalis' 01010101
sed -i '17i\>Calopteryx_splendens' 01010101
sed -i '19i\>Campodea_augens' 01010101
sed -i '21i\>Capnura_wanica' 01010101
sed -i '23i\>Catajapyx_aquilonaris' 01010101
sed -i '25i\>Cloeon_dipterum' 01010101
sed -i '27i\>Coccinella_septempunctata' 01010101
sed -i '29i\>Daphnia_magna' 01010101
sed -i '31i\>Galloisiana_yuasai' 01010101
sed -i '33i\>Inocellia_crassicornis' 01010101
sed -i '35i\>Laupala_kohalensis' 01010101
sed -i '37i\>Lepidocampa_weberi' 01010101
sed -i '39i\>Lipothrix_lubbocki' 01010101
sed -i '41i\>Machilis_hrabei' 01010101
sed -i '43i\>Mantispa_styriaca' 01010101
sed -i '45i\>Meinertellus_cundinamarcensis' 01010101
sed -i '47i\>Neelides_sp' 01010101
sed -i '49i\>Occasjapyx_japonicus' 01010101
sed -i '51i\>Octostigma_sinensis' 01010101
sed -i '53i\>Oncopodura_yosiiana' 01010101
sed -i '55i\>Oropsylla_silantiewi' 01010101
sed -i '57i\>Pediculus_humanus' 01010101
sed -i '59i\>Penaeus_vannamei' 01010101
sed -i '61i\>Plectrocnemia_conspersa' 01010101
sed -i '63i\>Protohermes_xanthodes' 01010101
sed -i '65i\>Pseudachorutes_palmiensis' 01010101
sed -i '67i\>Pseudobourletiella_spinata' 01010101
sed -i '69i\>Ptilocerembia_catherinae' 01010101
sed -i '71i\>Rhopalosiphum_maidis' 01010101
sed -i '73i\>Sinella_curviseta' 01010101
sed -i '75i\>Sinentomon_erythranum' 01010101
sed -i '77i\>Stylops_melittae' 01010101
sed -i '79i\>Tanzaniophasma_sp' 01010101
sed -i '81i\>Thalassaphorura_encarpata' 01010101
sed -i '83i\>Thermobia_domestica' 01010101
sed -i '85i\>Thrips_palmi' 01010101
sed -i '87i\>Tigriopus_japonicus' 01010101
sed -i '89i\>Timema_cristinae' 01010101
sed -i '91i\>Tricholepidion_gertschi' 01010101
sed -i '93i\>Zootermopsis_nevadensis' 01010101
sed -i '95i\>Zorotypus_caudelli' 01010101

#最后手动修改文件名即可