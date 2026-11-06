##### BioGeoBEARS: estimativa de áreas ancestrais (DEC, DIVALIKE, BAYAREALIKE, cada um com e sem +J) #####
# Versão traduzida e comentada em português do script biogeobears.R (a versão em inglês é a de referência)
# Adaptado dos scripts de exemplo do BioGeoBEARS de N. J. Matzke
# (http://phylo.wikidot.com/biogeobears; Matzke 2013, Frontiers of Biogeography 5:242-248; Matzke 2014, Systematic Biology 63:951-970)
# Instale os pacotes com install_packages.R e rode este script a partir da pasta biogeography/. Licença: GPL-3

dir.create("results_pt", showWarnings=FALSE)

##### Script BioGeoBEARS #####

# Para instalar o BioGeoBEARS é necessário instalar o Rtool
# (a forma de instalar o Rtool pode mudar conforme a versão do R)

# Instalar:


# Abrir:

library(Rcpp)
library(RcppArmadillo)
library(gdata)
library(gtools)
library(xtable)
library(plotrix)
library(vegan)
library(FD)
library(SparseM)
library(ape)
library(rexpokit)
library(cladoRcpp)
library(GenSA)
library(hms)
library(phytools)
library(expm)
library(fdrtool)
library(statmod)
library(spam)
library(MultinomialCI)
library(optimx)
library(snow)
library(phylobase)
library(devtools)

#A instalação do BioGeoBEARS é feita utilizando o devtools e É direto do github:


#  O pacote foi desenvolvido para inferências de ML e Bayesian de
# (a) distâncias geográficas ancestrais e
# (b) modelo para a evolução do alcance geográfico ao longo da filogenia

# Abrir:

library(BioGeoBEARS)

# Para definir o diretório principal:

# Run this script from the biogeography/ folder

# Os arquivos extras dos pacotes ficam no diretório "extdata" dentro dos pacotes instalados
# BioGeoBEARS contém muitos arquivos de exemplos no diretório extdata
# Para encontrar o diretório extdata

extdata_dir = normalizePath(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# A função "system.file" procura o diretório dentro do pacote indicado
# nesse caso o pacote indicado foi o BioGeoBEARS
# A função plot_BioGeoBEARS_results exige um script do diretório extdata
# para calcular as posições dos "corners" na plotagem.
# Essa função não roda direto no BioGeoBEARS porque usa rotinas do C
# do pacote APE que não passa na checagem do R CMD.

# Definir os arquivos da filogenia e distribuição:

trfn = "input/tree.nwk"
geogfn = "input/geography.data"

# Arquivo da filogenia #
#  1.  Precisa ser binário/bifurcando: sem politomias.
#  2.  Sem branchlengths negativos(e.g. BEAST MCC árvores consenso às vezes têm branchlengths negativos).
#  3.  Cuidado com ramos curtos, BioGeoBEARS irá interpretar ramos muito curtos como # ancestrais diretos.
#  4. Você pode usar árvores não-ultrametricas, mas BioGeoBEARS irá interpretar qualquer ponto significativamente abaixo do topo da árvore como fóssil.
#  Isso só é uma boa ideia se você possui ancestrais na sua árvore. Ver e.g. Wood, Matzke et al. (2013), Systematic Biology.
#  5.  O padrão do BioGeoBEARS faz sentido para árvores com branchlengths em unidades de milhões de anos e a árvore é 1-1000 units tall.
# Se tiver uma árvore com total height de e.g. 0.00001, você precisa ajustar e.g. os valores máximos d, e, ou (simpler) multiplique todos seus branchlengths
# para deixá-los em unidades razoáveis.
#  6. Não use espaços nos nomes das espécies, use E.G. "_"

# Verificar o arquivo raw Newick:

moref(trfn)

# Verifique a sua filogenia(plote em um PDF, para evitar muitos gráficos na mesma janela)

pdffn = "results_pt/tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(trfn)
tr
plot(tr)
title("ISIRA O TITULO")
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

dev.off()

# Arquivo Geográfico #
#  1. Arquivo formatado tipo PHYLIP.
# O que significa que na primeira linha
#    -  1st número é o número de linhas (especies/demes)
#    -  2nd número é o número de colunas (numéro de áreas)
#    -  depois do tab, coloque as áreas em parenteses, com espaços: (A B C D)
#
#  1.5. Exemplo primeira linha:
#    10    4    (A B C D)
#
#  2. A segunda linha e as seguintes:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
#
#  2.5a. TAB entre os nomes das espécies e a área 0/1s
#  2.5b. Sem espaço e sem TAB entre a área 0/1s.
#  3. Todos os nomes no arquivo de geografia precisam ser iguais aos nomes no arquivo de filogenia.
#  4. Não use espaços nos nomes das espécies, use E.G. "_"
#  5. Unidades taxonomicas operacionais (OTUs) deveriam ser linhagens filogenéticas i.e. populações geneticamente isoladas.
#Podendo ou não ser idênticos com as espécies. Você não quer usar apenas espécimes, como cada espécime
#pode viver apenas em uma área, o que vai favorecer os modelos DEC+J
#Não tem problema se as espécies/linhagens realmente vivem em uma única área mas você não vai querer assumir isso sem pensar muito a respeito.
#Em resumo, você deve collapsar multiplos especimens em espécie/linhagem se dados indicam que são a mesma população genética.

# Verificar o arquivo raw geography:
moref(geogfn)

# Verificar a abrangência dos dados geográficos:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Tamanho máximo do alcance observado:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Indique o máximo de áreas que qualquer espécie pode ocupar
# não pode ser maior que o número de áreas que você indicar, mas pode ser menor.

max_range_size = 12

# DICA: O número de estados (= número de ranges geográficos possíveis)
# depende do (a) número de áreas e o (b) max_range_size.
# Se você tem mais de 500-600 estados, os cálculos serão muito lentos
# pois o programa exponencia a matriz de e.g. 600x600.
# O computador pode travar e nem completar a primeira likelihood.
# Geralmente é o que ocorre quando o LAGRANGE diz: you have too many states!
# Para chegar o número de estados para determinada abrangência:
#numstates_from_numareas(numareas=12, maxareas=9, include_null_range=TRUE)

# Um grande número de áreas gera problemas:...exceto se você limitar o max_range_size:
#numstates_from_numareas(numareas=12, maxareas=4, include_null_range=TRUE)

#################################################################################
# DEC e DEC+J #
#################################################################################
# O modelo DEC do BioGeoBEARS é idêntico ao modelo do LAGRANGE
# e deve retornar valores idênticos de ML, log-likelihoods e datasets.
#
# Probabilidades de estados ancestrais dos nós serão levemente diferentes
# uma vez que o BioGeoBEARS reporta as probabilidades de estados ancestrais
# em um modelo ML global e o LAGRANGE reporta após re-otimizar a likelihood
# após ajeitar/fixar o estado de cada nó. Podem ser similares, mas não idênticos.
# Matzke (2014), Systematic Biology

# Inicie o modelo padrão (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Forneça ao BioGeoBEARS a localicação do arquivo phylogeny Newick
BioGeoBEARS_run_object$trfn = trfn

# Forneça ao BioGeoBEARS a localicação do arquivo geography text
BioGeoBEARS_run_object$geogfn = geogfn

# Informe o maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min to treat tip as a direct ancestor (no speciation event)
# BioGeoBEARS_run_object$min_branchlength = 0.000001

# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$include_null_range = TRUE

# Definindo uma análise de tempo estratificada (arquivos em input/):
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Se necessário use mais núcleos para acelerar as análises, mas isso requer o
# library(parallel) e/ou library(snow).
# Opções de velocidade e processamento multicore
BioGeoBEARS_run_object$on_NaN_error = -1e50    # retorna um valor baixo de lnL se os parametros produzirem um erro NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # atalhos para acelerar a busca ML; use FALSE se preocupar (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # se FALSE, use optim() no lugar de optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1

# Esta função carrega a matriz do multiplicador de dispersão etc. dos arquivos de texto # no objeto de modelo. É necessário para que funcionem!
# (Ele também executa algumas verificações nos inputs para determinados erros.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide a árvore por timeperiods/strata (uncomment para análise estratificada)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# A árvore estratificada é descrita nessa tabela:
BioGeoBEARS_run_object$master_table

# Boas configurações para obter estados ancestrais:
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # obter os estados ancestrais da corrida optim

# Definir o modelo DEC (nada a fazer; defaults)

# Veja o BioGeoBEARS_run_object; é só uma lista de configurações.
BioGeoBEARS_run_object

# Veja o modelo dos objetos
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# Veja a tabela que contém os parâmetros do modelo
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Rode esse para confirmar os inputs. Considere as mensagens de erro!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Para uma análise lenta, corra uma vez, então defina runslow=FALSE para rodar apenas o resultado salvo.

runslow = TRUE
resfn = "results_pt/DEC_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)
    resDEC = res
    } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
    }

#################################################################################################################################
# Rodando DEC+J
#################################################################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min para considerar a tip como ancestral direto (sem evento de especiação)
BioGeoBEARS_run_object$include_null_range = TRUE    # FALSE para e.g. DEC* model, DEC*+J, etc.

# Defina a análise de tempo estratificado
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Se necessário use mais cores para acelerar as analises
BioGeoBEARS_run_object$on_NaN_error = -1e50    # retorna um baixo lnL se os parametros produzirem erro NaN (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # atalhos para acelerar a busca ML; use FALSE se preocupar (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # se FALSE, use optim() no lugar de optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causa patologias e não é muito mais rápido
BioGeoBEARS_run_object$num_cores_to_use = 1

# Essa função carrega a matrix do multiplicador de dispersão etc. do arquivo de texto no objeto modelo.
# É necessária para funcionar! (Ela também executa algumas verificações nos inputs para determinados erros.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide a árvore por timeperiods/strata (uncomment para análise estratificada)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# A matrix estratificada está descrita na tabela:
BioGeoBEARS_run_object$master_table

# Boas configurações para conseguir estados ancestrais
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Definir o modelo DEC+J
# Consegue os valores do parametro ML do modelo nested 2-parametros
# (isso garante que o modelo de 3-parametros sempre será pelo menos bom)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Valores iniciais do Input para d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Adicione o j como parâmetro livre
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "results_pt/DEC+J_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)

    resDECj = res
    } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
    }


# PDF plots
pdffn = "results_pt/DEC+J_v1.pdf"
pdf(pdffn, width=6, height=6)

###############################################################################################
# Plotando o estado ancestral - DEC #
###############################################################################################
analysis_titletxt ="COLOCAR O NOME QUE DESEJA"

# Configuração
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Estados
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Gráfico de pizza
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

###############################################################################################
# Plotando estados ancestrais - DECJ #
###############################################################################################
analysis_titletxt ="COLOCAR O NOME QUE DESEJA"

# Configuração
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Estados
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Gráfico de pizza
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF

##############################################################################################################
# DIVALIKE e DIVALIKE+J #
##############################################################################################################
# O modelo "DIVALIKE" do BioGeoBEARS não é identico ao modelo de Ronquist (1997), DIVA com parcimonia.
# É uma interpretação de verosimilhança do DIVA, construindo os processos do modelo DIVA do mesmo modo
# que o modelo DEC, mas somente permitindo processos que o DIVE permite
# (vicariancia difundida: sim; subset de simpatria: não; veja Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE é uma interpretação de verossimilhanã da DIVA de parcimonia, é parecido com a DIVA
# similar, não identico, por isso se chama "DIVALIKE" ;-)
##############################################################################################################

##############################################################################################################
# Rodando DIVALIKE #
##############################################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min para considerar a tip como um ancestral direto (sem evento de especiação)
BioGeoBEARS_run_object$include_null_range = TRUE    # FALSE para e.g. DEC* model, DEC*+J, etc.

# Definir as análises de tempo estratificado:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Opções para acelerar e processor multicore
BioGeoBEARS_run_object$on_NaN_error = -1e50    # retorna um baixo lnL se os parametros produzem um erro NaN (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # atalhos para acelerar a busca ML; use FALSE se preocupar (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # se FALSE, use optim() no lugar de optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causa patologia e não é tão mais rápido
BioGeoBEARS_run_object$num_cores_to_use = 1

# Essa função carrega a matrix multiplicadora de dispersão etc. do arquivo de texto para o modelo objeto.
# Isso é necessário para funcionar!(Também correm algumas checagens para certos erros desses inputs)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide a árvore pelo timeperiods/strata (uncomment para a análise estratificada)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# A arvore estratifica é descita na tabela:
BioGeoBEARS_run_object$master_table

# Bons padrões de configuração para os estados ancestrais
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # obtem os estados ancestrais da corrida optim

# Configurando o modelo DIVALIKE
# Sem a simpatrica com subset
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Permite a vicariancia classica, difundidade, todos os eventos com mesma probabilidade
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Sem salto de dispersão/evento fundador
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "results_pt/DIVALIKE_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)
    resDIVALIKE = res
    } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKE = res
    }

#######################################################
# Rodando DIVALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min para considerar a tip como um ancestral direto (sem evento de especiação)
BioGeoBEARS_run_object$include_null_range = TRUE    # FALSE para e.g. DEC* model, DEC*+J, etc.

# Configurando a análise de tempo estratificado:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Opções para acelerar e processor multicore
BioGeoBEARS_run_object$on_NaN_error = -1e50    # retorna um baixo lnL se os parametros produzem um erro NaN (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # atalhos para acelerar a busca ML; use FALSE se preocupar (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # se FALSE, use optim() no lugar de optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causa patologia e não é tão mais rápido
BioGeoBEARS_run_object$num_cores_to_use = 1

# Essa função carrega a matrix multiplicadora de dispersão etc. do arquivo de texto para o modelo objeto.
# Isso é necessário para funcionar!(Também correm algumas checagens para certos erros desses inputs)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide a árvore pelo timeperiods/strata (uncomment para a análise estratificada)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# A arvore estratificada é descita na tabela:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # obtem estados ancestrais para a corrida optim

# Configura o modelo DIVALIKE+J
# Obtem os valores de parametros ML do modelo nested de 2-parametros
# (isso garante que o modelo de 3-parametros seja sempre ao menos bom)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Valores iniciais para d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Sem simpatrica com subset
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Permite a vicariancia classica, difundida, com todos eventos com mesma probabilidade
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Adiciona o salto de dispersão/Evento fundador
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# No DIVALIKE+J, o máximo de "j" deve ser 2, não 3 (como é o padrão do DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "results_pt/DIVALIKE+J_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)

    resDIVALIKEj = res
    } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKEj = res
    }

pdffn = "results_pt/DIVALIKE_vs_DIVALIKE+J_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plotando os estados ancestrais - DIVALIKE
#######################################################
analysis_titletxt ="COLOCAR O NOME QUE DESEJA"

# Configurando
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Estados
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Gráfico de pizza
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plotando os estados ancestrais - DIVALIKE+J
#######################################################
analysis_titletxt ="COLOCAR O NOME QUE DESEJA"

# Configurando
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Estados
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Gráfico de pizza
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

###############################################################################################
# BAYAREALIKE e BAYAREALIKE+J #
###############################################################################################
# Assim como no DIVA, o modelo BioGeoBEARS BayArea-like não é
# idêntico com o modelo Bayesiano implementado por completo no
# programa "BayArea" de Landis et al. (2013).
#
# Nesse caso é uma interpretação likelihood simplificada do modelo
# Basicamente na BayArea e no BioGeoBEARS-BAYAREALIKE, "d" e "e" funcionam como no modelo
# DEC do Lagrange (e BioGeoBEARS), e então a cladogeneses assumida pela BayArea's
# não é nada particular o que acontece na cladogeneses) é replicado pelo BioGeoBEARS.
#
# O que nos deixam 3 coisas importantes na BayArea:
# 1. Dependência da distância (você pode adicionar isso com uma
# matrix de distância + o parâmetro "x" no BioGeoBEARS)
# 2. Uma correção para desautoriza os eventos "e" que levam a espécie à extinção
# (um raio geográfico nulo)
# 3. A bela amostragem bayesiana de histórias, que permite análises em um grande
# número de áreas.
# O foco principal de um modelo "BAYAREALIKE" é para testar a importância do
# modelo de cladogeneses em bancos de dados específicos.  Ajuda ou prejudica
# os dados de likelihood se não existir um processo especial de cladogeneses?
#
# BAYAREALIKE é uma interpretação de likelihood da BayArea, e é "tipo BayArea"
# Similar, mas não idêntico ao Bayesian BayArea. Por isso chamo o modelo de "BAYAREALIKE".
###############################################################################################

###############################################################################################
# Rodando o BAYAREALIKE
###############################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min para considrar a tip como ancestral direto (sem evento de especiação)
BioGeoBEARS_run_object$include_null_range = TRUE    # FALSE para e.g. DEC* model, DEC*+J, etc.

# Configurando uma análise de tempo estratificado:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Opções para acelerar e processor multicore
BioGeoBEARS_run_object$on_NaN_error = -1e50    # retorna um baixo lnL se os parametros produzem um erro NaN (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # atalhos para acelerar a busca ML; use FALSE se preocupar (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # se FALSE, use optim() no lugar de optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causa patologia e não é tão mais rápido
BioGeoBEARS_run_object$num_cores_to_use = 1

# Essa função carrega a matrix multiplicadora de dispersão etc. do arquivo de texto para o modelo objeto.
# Isso é necessário para funcionar!(Também correm algumas checagens para certos erros desses inputs)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide a árvore pelo timeperiods/strata (uncomment para a análise estratificada)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# A árvore estratificada pode ser descrita nessa tabela:
BioGeoBEARS_run_object$master_table

# Bons padrões de configurações para conseguir estados ancestrais
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # obtem os estados ancestrais para a corrida optim

# Configure o modelo BAYAREALIKE
# Sem subset de simpatria
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# Sem vicariancia
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Sem especiação por dispersão/efeito-fundador
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Ajuste a ligação entre os paramêtros
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Somente eventos simpatricos/range-copying (y) permitidos, e com cópia exata
# (ambos descendentes sempre dos mesmo tamanho que o ancestral)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Conferindo os inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "results_pt/BAYAREALIKE_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)
    resBAYAREALIKE = res
    } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKE = res
    }

###################################################################################################################
# Rodando BAYAREALIKE+J
###################################################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min para considrar a tip como ancestral direto (sem evento de especiação)
BioGeoBEARS_run_object$include_null_range = TRUE    # FALSE para e.g. DEC* model, DEC*+J, etc.

# Configurando uma análise de tempo estratificado:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Opções para acelerar e processor multicore
BioGeoBEARS_run_object$on_NaN_error = -1e50    # retorna um baixo lnL se os parametros produzem um erro NaN (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # atalhos para acelerar a busca ML; use FALSE se preocupar (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # se FALSE, use optim() no lugar de optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causa patologia e não é tão mais rápido
BioGeoBEARS_run_object$num_cores_to_use = 1

# Essa função carrega a matrix multiplicadora de dispersão etc. do arquivo de texto para o modelo objeto.
# Isso é necessário para funcionar!(Também correm algumas checagens para certos erros desses inputs)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide a árvore pelo timeperiods/strata (uncomment para a análise estratificada)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# A árvore estratificada pode ser descrita nessa tabela:
BioGeoBEARS_run_object$master_table

# Bons padrões de configurações para obter os estados ancestrais
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # obtem estados ancestrais da corrida optim

# Configurando o modelo BAYAREALIKE+J
# Conseguindo os valores do parametro ML do modelo nested de 2-parametros
# (isso certificará que o modelo 3-parameter sempre será tão bom quanto)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Configurando os valores iniciais de d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Sem subset de simpatria
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# Sem vicariancia
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Permitir o jump dispersal/founder-event speciation
# (configurando o valor inicial próximo de 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# No BAYAREALIKE+J, o máximo de "j" deve ser 1
# não 3 (como no padrão de DEC+J)ou 2 (como no DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Ajustando a linkage entre os parametros
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Somente eventos sympatric/range-copying (y) permitidos, e com
# cópia exata (ambos descendentes sempre com mesmo tamanho que o ancestral)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# BAYAREALIKE+J ás vezes trava em alguns computadores, principalmente Windows
# Algumas precisões no problema da corrida, quando optim/optimx tenta alguns valores
# de parâmetros abaixo de 0. O “min” e “max” de cada parâmetro deveriam prevenir
# isso, mas aparentemente o optim/optimx ultrapassa esses limites. Caso o computador travar,
# tente aumentar um pouco o "min" e diminuir um pouco o "max" de cada parâmetro
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "results_pt/BAYAREALIKE+J_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)

    resBAYAREALIKEj = res
    } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKEj = res
    }

pdffn = "results_pt/BAYAREALIKE_vs_BAYAREALIKE+J_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plotando os estados ancestrais - BAYAREALIKE
#######################################################
analysis_titletxt ="COLOCAR O NOME DESEJADO"

# Configuração
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Estados
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Gráfico de pizza
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plotando os estados ancestrais - BAYAREALIKE+J
#######################################################
analysis_titletxt ="COLOCAR O NOME DESEJADO"

# Configurando
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Estados
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Gráfico de pizza
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

#######################################################################
# Configurando as tabelas vazias para incluir os resultados estatísticos
restable = NULL
teststable = NULL

#######################################################################
# Estatística - DEC vs. DEC+J
#######################################################################
# Nós temos que extrair o log-likelihood de forma diferente. Dependendo da versão do optim/optimx.
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, modelo nulo para o Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# A hipotese nula para o Likelihood Ratio Test (LRT) é de que 2 modelos conferem o mesmo likelihood nos dados.

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Estatistica -- DIVALIKE vs. DIVALIKE+J
#######################################################
# Nós temos que extrair o log-likelihood de forma diferente. Dependendo da versão do optim/optimx.
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, modelo nulo para o Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# DIVALIKE+J, modelo alternativo para o Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Estatistica -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# Nós temos que extrair os log-likelihood de formas diferentes
# Dependendo da versão of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, modelo nulo para o Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################################
# Reunindo as tabelas de resultados: DEC, DEC+J, DIVALIKE, DIVALIKE+J, ...
#######################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# Veja os resultados!!
restable
teststable

#######################################################
# Salve os resultados das tabelas para depois – cheque para e.g. problemas convergentes
#######################################################

# Carregando na "restable"
save(restable, file="results_pt/restable_v1.Rdata")
load(file="results_pt/restable_v1.Rdata")

# Carregando na "teststable"
save(teststable, file="results_pt/teststable_v1.Rdata")
load(file="results_pt/teststable_v1.Rdata")

# Também salve nos arquivos de texto
write.table(restable, file="results_pt/restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="results_pt/teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Peso do Modelo para todos os seis modelos
#######################################################
restable2 = restable

# Com AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# Com AICcs – fatores no tamanho da amostra
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

# Também salve em arquivos de texto
write.table(restable_AIC_rellike, file="results_pt/restable_AIC_rellike.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="results_pt/restable_AICc_rellike.txt", quote=FALSE, sep="\t")

# Salve com boas condições de formatação
write.table(conditional_format_table(restable_AIC_rellike), file="results_pt/restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="results_pt/restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\t")

writeLines(capture.output(sessionInfo()), "results_pt/sessionInfo.txt")
