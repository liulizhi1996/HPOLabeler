# HPOLabeler

## 如何运行？

**注：请在开始下面步骤前，在~/.bashrc文件的最后一行加上：**

	export PYTHONPATH=${PYTHONPATH}:[HPOLabeler的目录位置]

---

### 第一步：预处理（Pre-processing）

1. 从[http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/](http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/)下载基因-HPO标注文件`ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt`到目录`data/annotation/raw`。然后，运行`src/preprocessing/extract_gene_id.py`以提取文件的第一列（即entrez gene id），并将其写入一个.txt文件。接着，上传该文件内容到Uniprot ID Mapping Tool [(http://www.uniprot.org/mapping/)](http://www.uniprot.org/mapping/)，选择选项`From: Entrez Gene (Gene ID)`和`To: UniProtKB`，然后点击`Submit`按钮。在新的窗口中，选择左栏的`Filter by`为`Reviewed Swiss-Prot`，并点击页面中间的`Download`按钮，选择`Format: Tab-separated`并点击`Go`下载映射文件。将映射文件放入`data/annotation/intermediate`中，并更名为`gene2uniprot.txt`。

2. 运行`src/preprocessing/create_annotation.py`程序。程序将输出处理好的HPO标注文件。

3. 重复上述两步，分别处理得到三个标注文件（比如2018-03-09、2019-04-15和2019-11-15）。

4. 从[https://bioportal.bioontology.org/ontologies/HP](https://bioportal.bioontology.org/ontologies/HP)下载三个时期相对应的.obo文件，将其放在`data/obo`目录下。之后运行`src/preprocessing/split_dataset.py`。我们将得到：
	- 处理好的用于训练基础分类器的标注、用于训练排序学习的标注和用于测试的标注
	- 上面三个标注文件中的蛋白质标识符列表
	- 用来标注蛋白质的HPO term列表
	- 按照频率划分的各个小组内的HPO term列表

### 第二步：处理特征（Extracting Features）

#### STRING

1. 打开[https://string-db.org/](https://string-db.org/)，然后点击页面左上角的`Version`，在自动跳转到的新页面中选择合适的版本，并单击`Address`一栏中的链接。之后，点击新页面上方导栏的`Download`按钮，点击`choose an organism`下拉菜单，选择`Homo sapiens`。现在，点击`INTERACTION DATA`部分的`9606.protein.links.XXX.txt.gz`（XXX为版本），下载蛋白质互作数据。最后，点击`ACCESSORY DATA`部分的`mapping_files (download directory)`，进入ftp页面，点击`uniprot_mappings/`目录，下载属于人类（可能是开头为`9606`或者文件名中有human字样）的压缩文件。上述两个文件都下载至`data/feature/STRING/raw`目录下。

2. 运行`src/feature/STRING/string.py`程序，得到在`data/feature/STRING/clean`目录下的json格式的文件。

#### GeneMANIA

1. 从[http://genemania.org/data/current/Homo\_sapiens.COMBINED/](http://genemania.org/data/current/Homo_sapiens.COMBINED/)下载`COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt`文件。然后再从[http://genemania.org/data/current/Homo_sapiens/](http://genemania.org/data/current/Homo_sapiens/)中下载`identifier_mappings.txt`文件。

2. 运行`src/feature/GeneMANIA/genemania.py`，得到在`data/feature/GeneMANIA/clean`下的json格式文件。

#### BioGRID

1. 从[https://downloads.thebiogrid.org/BioGRID/Release-Archive](https://downloads.thebiogrid.org/BioGRID/Release-Archive)中选取合适的版本，并点击链接进入。在新页面中下载`BIOGRID-ALL-XXX.tab2.zip`（XXX为版本号）。然后在[https://downloads.thebiogrid.org/BioGRID/External-Database-Builds/](https://downloads.thebiogrid.org/BioGRID/External-Database-Builds/)中下载`UNIPROT.tab.txt`。这两个文件都放到`data/feature/BioGRID/raw`下。

2. 运行`src/feature/BioGRID/biogrid.py`，得到在`data/feature/BioGRID/clean`下的json文件。

#### GO Annotation

1. 首先，在[https://www.ebi.ac.uk/GOA/downloads](https://www.ebi.ac.uk/GOA/downloads)页面的`Annotation Sets`表格的`Human`一行中，点击某一个链接（若要下载当前最新数据，点击`Current Files`，否则点击`Archive Files`），然后下载合适版本的.gaf文件。

2. 若要下载最新版的GO的.obo文件，可从[http://geneontology.org/docs/download-ontology/](http://geneontology.org/docs/download-ontology/)中下载；若要下旧版，则可以从[https://bioportal.bioontology.org/ontologies/GO](https://bioportal.bioontology.org/ontologies/GO)下载。

3. 运行`src/feature/GO_annotation/go_annotation.py`，得到在`data/feature/GO_annotation/clean`目录下的数据。

#### Trigram

1. 将`data/dataset/protein`中的三个蛋白质标识符列表中所有的蛋白质标识符上传至[https://www.uniprot.org/mapping/](https://www.uniprot.org/mapping/)，并且选择`From`为`UniProtKB AC/ID`，`To`为`UniProtKB`，再点击`Submit`。在新页面中，点击页面中间的`Download`按钮，选择格式为`FASTA (canonical)`，点击`Go`。将下载得到的fasta序列文件放到`data/feature/Trigram/raw`中。

2. 运行`src/feature/Trigram/trigram.py`，得到在`data/feature/Trigram/clean`中的输出文件。

#### InterPro

1. 从[http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/)选择合适的InterProScan程序包下载。

2. 进入下载后解压的目录内，以要查询的蛋白质的fasta文件为输入，运行InterProScan，得到程序匹配到的InterPro signatures的xml文件：
	
	`./interproscan.sh -i /path/to/sequences.fasta -b /path/to/output_file -f XML`

3. 运行`src/feature/InterPro/interpro.py`程序，处理上一步得到的原始xml文件，获得处理后的InterPro特征文件。

#### HumanNet

1. 打开[https://www.inetbio.org/humannet/download.php](https://www.inetbio.org/humannet/download.php)。选择`Integrated Networks`下的`HumanNet-XN`下载。将下载到的文件放在`data/feature/HumanNet/raw`中。

2. 打开[https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22](https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22)。点击页面中部的`Columns`按钮，在新页面的`Columns to be displayed`中点击所有的虚线框右上角的叉号。之后在`Add more columns`栏目中部的`Search:`搜索框中输入`GeneID`，点击弹出的联想词。此时，单击最右侧的`Save`按钮。跳回到原来的页面，这时表格只剩下`Entry`和`Cross-reference (GeneID)`。点击页面中间的`Download`按钮，选择`Format: Tab-separated`，再点击`Go`按钮下载文件。将该文件重命名为`entrez2uniprot.txt`，放在`data/feature/COXPRESdb/raw`下。下载完毕后，在刚刚UniProt的页面中，再次点击`Columns`按钮，在新页面中点击右侧的`Reset to default`（注意：不要点击default这个单词，而点击Reset），之后再单击`Save`。这样UniProt界面恢复原来的样子。

3. 运行`src/feature/HumanNet/humannet.py`，得到处理后的蛋白质关系网络。

#### HIPPIE

1. 打开[http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)，并选择HIPPIE tab format下的合适版本进行下载。

2. 打开[https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22](https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22)。点击页面中部的`Columns`按钮，在新页面的`Columns to be displayed`中点击除`Entry name`外所有的虚线框右上角的叉号，再单击最右侧的`Save`按钮。跳回到原来的页面，此时表格只剩下`Entry`和`Entry name`两列。点击页面中间的`Download`按钮，选择`Format: Tab-separated`，再点击`Go`按钮下载文件。将该文件重命名为`uniprot_name2id.txt`，放在`data/feature/HIPPIE/raw`下。下载完毕后，在刚刚UniProt的页面中，再次点击`Columns`按钮，在新页面中点击右侧的`Reset to default`（注意：不要点击default这个单词，而点击Reset），之后再单击`Save`。这样UniProt界面恢复原来的样子。

3. 运行`src/feature/HIPPIE/hippie.py`，获得处理后的HIPPIE的PPI数据。

#### COXPRESdb

1. 打开[https://coxpresdb.jp/download/](https://coxpresdb.jp/download/)，点击`> Show previous data`。找到`Species`列为`Human`，`Method`列为`R`的行，下载这一行`Gene correlation table`列的.zip文件，将其放入`data/feature/COXPRESdb/raw`目录下。**注意：请不要解压！**

2. 打开[https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22](https://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22)。点击页面中部的`Columns`按钮，在新页面的`Columns to be displayed`中点击所有的虚线框右上角的叉号。之后在`Add more columns`栏目中部的`Search:`搜索框中输入`GeneID`，点击弹出的联想词。此时，单击最右侧的`Save`按钮。跳回到原来的页面，这时表格只剩下`Entry`和`Cross-reference (GeneID)`。点击页面中间的`Download`按钮，选择`Format: Tab-separated`，再点击`Go`按钮下载文件。将该文件重命名为`entrez2uniprot.txt`，放在`data/feature/COXPRESdb/raw`下。下载完毕后，在刚刚UniProt的页面中，再次点击`Columns`按钮，在新页面中点击右侧的`Reset to default`（注意：不要点击default这个单词，而点击Reset），之后再单击`Save`。这样UniProt界面恢复原来的样子。

3. 运行`src/feature/COXPRESdb/coxpresdb.py`，获得处理后的共表达数据。

### 第三步：训练基础分类器（Basic Models）

#### Naive

1. 运行`src/basic/naive/naive.py`，得到`data/result/basic/naive`内的输出文件。

#### Neighbor

1. 注意设置配置文件里“network”一项的“type”，对于STRING和GeneMANIA，其设为“weighted”；对于BioGRID，则要设置为“unweighted”。

2. 运行`src/basic/neighbor/neighbor.py`，得到保存在`data/result/basic/neighbor`中的预测结果。

#### Flat

1. 运行`src/basic/flat/flat.py`，将各种处理得到的特征文件作为输入，训练Logistic Regression分类器，对用于排序学习的训练集和测试集进行预测，得到输出在`data/result/basic/flat`目录下的一系列预测结果文件。

### 第四步：排序学习（Learning to Rank）

1. 将第三步中得到的一系列预测分数作为输入，即配置文件中的`"result"`部分。**注意：请务必保证这一部分的`"ltr"`和`"test"`的列表内的文件顺序是一致的！**

2. 注意合理调节配置文件中`"model_param"`部分的`max_depth`值、`"fit_param"`部分中的`num_boost_round`和`early_stopping_rounds`以及`"top"`中的取值。

3. 运行`src/ensemble/ltr/ltr.py`程序，得到最终的预测结果。


### 第五步：评估（Evaluation）

1. 将要评估的预测结果的文件路径添加在配置文件的`"result"`部分。

2. 运行`src/utils/evaluation.py`，程序将会显示各个预测结果在各个子本体上的
	- Fmax：以蛋白质为中心的评估指标
	- AUC：以HPO term为中心的评估指标，即每个HPO term的AUC的平均值
	- AUPR：整体的评估指标，即以一对蛋白质-HPO term为实例进行计算的AUPR
	- 每个按频率划分的HPO term小组内的平均AUC

