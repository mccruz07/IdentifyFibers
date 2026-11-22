# IdentifyFibers RC2 · 2D / 3D

Candidato privado para testes no CellProfiler 4.2.8. O mesmo módulo oferece **2D image** e **3D volume (test candidate)**. Um pipeline salvo na v1.0.1 abre em 2D. O modo 3D usa realce tubular Sato, extração de linhas centrais e continuação direcional inspirada no FIRE. A extração volumétrica tem lógica própria.

[English guide](README.md) · [Cinco auditorias e resultados](AUDIT.md) · [Manual do modo 2D](docs/manual.pt-BR.md)

## Como testar a atualização 2D

Escolha **2D image** e abra o mesmo pipeline com a mesma imagem e os mesmos parâmetros nas duas versões. Cronometre o IdentifyFibers com a janela de visualização do módulo fechada. Compare os rótulos, as linhas centrais e as medidas exportadas. [Tempos e validação](PERFORMANCE.md)

## Primeiro teste

1. Crie uma pasta de plugins para o teste. Copie para ela o arquivo `identifyfibers.py` e a pasta completa `identifyfibers_engine`, ambos em `active_plugins/`. Selecione essa pasta nas preferências do CellProfiler e reinicie. Guarde o caminho da instalação habitual para poder voltar à v1.0.1.
2. Abra **demo/Fiber3D_demo.cppipe** e adicione **apenas demo/fiber_phantom.ome.tif** à lista de imagens. Escolha uma pasta de saída. Em **Test Mode**, percorra os planos Z na janela do IdentifyFibers. Use **Analyze Images** para gravar os CSVs.
3. O volume sintético contém três fibras. Duas se cruzam na projeção XY e passam em planos diferentes. O resultado esperado é de três objetos, com comprimentos aproximados de **19,25, 20,50 e 23,69 µm**. O arquivo `demo/expected_fibers.csv` contém as medidas obtidas no teste completo.

![Projeção e profundidade do volume de teste](demo/phantom_views.png)

## Primeira pilha confocal

Use um único canal de fibras claras, um único tempo e um volume em Z. Em **NamesAndTypes**, marque **Process as 3D** e informe a razão entre os tamanhos dos voxels. Exemplo: com pixels de 0,25 µm em X e Y e passo de 0,80 µm em Z, use **X = 1, Y = 1, Z = 3,2**. Para **Y pixel size (µm)** no IdentifyFibers, insira **0,25**. O valor inicial zero pede que a calibração seja informada. Use o espaçamento entre os planos adquiridos, em micrômetros. Confira no CSV da imagem os três valores `Fiber3D_VoxelSize…_um`. O pipeline de demonstração já contém a calibração do exemplo. Ajuste-a ao trocar de aquisição.

| Controle | Como ajustar |
|---|---|
| Minimum / maximum expected width | Meça a espessura aparente de algumas fibras no sinal original. A faixa em µm seleciona três escalas de realce. Reduzir o mínimo favorece fibrilas finas e pode destacar ruído. Aumentar o máximo favorece estruturas mais largas e aumenta a distância usada para estimar direção nas junções. |
| Threshold correction factor | Comece em 1. Reduza se um trecho fraco e visível desaparece. Aumente se regiões de fundo entram no objeto. Examine vários planos para saber se o ajuste também une vizinhas. |
| Minimum traced length | Comprimento mínimo ao longo da curva em µm. Aumentar elimina caminhos curtos, inclusive ramos biológicos. |
| Maximum crossing deflection | Permite o desvio angular ao continuar de um ramo para outro. Valores baixos favorecem continuação reta. Aumentar aceita curvas mais fechadas e eleva a possibilidade de trocar de fibra. |
| CellProfiler binary volume | Use uma máscara binária registrada com a mesma forma e o mesmo espaçamento. Os voxels da máscara definem a segmentação. A entrada original fornece as intensidades. |
| Within objects | Cada objeto parental delimita um volume de análise. Uma fibra que cruza duas regiões gera trajetórias separadas. |

## Do volume às medidas

1. **Sinal e escala.** O algoritmo normaliza os voxels válidos por região e calcula uma resposta tubular com a matriz de segundas derivadas da intensidade (Hessiana), levando em conta a escala física nos três eixos. Duas curvaturas transversais negativas favorecem um centro tubular claro. Otsu aplicado às respostas positivas define a máscara. Quando o ruído é estimável, um teste adicional de contraste verifica o suporte do sinal original ao longo de cada trajetória. A máscara externa pula essas etapas.
2. **Trajetória.** Componentes usam 26 vizinhos. Linhas já finas seguem diretamente para a análise do grafo. Em estruturas espessas, o afinamento 3D de Lee fornece ramos e junções. Em cada junção, o algoritmo associa primeiro as direções fisicamente mais próximas de uma continuação reta. Cada ramo é usado uma única vez. Se o afinamento apagar um tubo, um caminho geodésico que privilegia o interior da máscara pode recuperar seu eixo. Essa alternativa exige alongamento global de ao menos 2:1, estimado pelos desvios nas direções principais; consulte `FallbackTrace`.
3. **Comprimento.** A trajetória é reamostrada por distância física e suavizada na escala de um passo entre voxels para reduzir o efeito de escada. Curvas abertas preservam as extremidades. Com `PathLength_um` obtém-se o comprimento dessa curva. `VoxelPathLength_um` fornece a soma bruta dos passos antes da suavização. Curvas muito fechadas e pontas pouco resolvidas requerem inspeção.
4. **Objetos.** As linhas aceitas servem de sementes. Dentro de cada componente com sementes, a distância percorrida em µm e as variações de intensidade determinam a atribuição de cada voxel a um objeto. O componente inteiro é distribuído; um prolongamento curto rejeitado ainda pode contribuir para o volume do objeto vizinho. Componentes sem trajetória aceita ficam como fundo.

## Como ler o CSV

As colunas abaixo começam com `Fiber3D_`. `EndToEnd_um` mede a distância reta entre as extremidades. `Straightness` divide essa distância por `PathLength_um`. `Volume_um3` soma o volume dos voxels atribuídos. `ApparentWidth_um` é a média de duas vezes a distância física do eixo ao fundo da máscara; a primeira camada fora do campo conta como fundo. Feixes e junções podem alargar a estimativa. `IntensityMean` e `IntensityIntegrated` usam os valores carregados pelo CellProfiler. A intensidade integrada é a soma dos valores dos voxels atribuídos, sem multiplicação pelo volume do voxel.

A direção principal da trajetória define `Azimuth_deg` no plano XY (X para a direita, Y no sentido crescente das linhas, módulo 180°) e `Inclination_deg`, inclinação assinada em relação a XY. Uma fibra vertical tem inclinação de 90° e azimute indefinido. Um eixo principal indeterminado produz ângulos ausentes. A orientação descreve a geometria, sem atribuir polaridade biológica. `Location_Center_X/Y/Z` indica o centroide da região em índices de voxel.

`TouchesImageEdge = 1` identifica uma região que alcança a borda do volume. `FallbackTrace = 1` identifica uma trajetória recuperada pelo caminho alternativo. Para investigar uma fibra, compare a imagem original, a máscara ou o realce e os rótulos ao percorrer Z. Um voxel tem um único dono, mas a identidade de fibras no interior de um feixe saturado pode permanecer indeterminada. A janela e as imagens de saída permitem inspecionar planos Z. A FFT pertence ao modo 2D.

**Escopo do modo 3D.** Os ensaios volumétricos desta distribuição usam dados sintéticos. O teste com microscopia incluído na suíte histórica é 2D. A validação biológica em pilhas reais permanece em aberto. O primeiro teste em dados próprios deve incluir uma fibra isolada, um feixe denso e um cruzamento aparente em XY. A distância entre planos e a resolução óptica axial limitam o que pode ser distinguido em profundidade. Comece com um recorte representativo e conserve o arquivo original, a calibração e o pipeline.

**Referência.** Stein AM, Vader DA, Jawerth LM, Weitz DA, Sander LM (2008). An algorithm for extracting the network geometry of three-dimensional collagen gels. *Journal of Microscopy*, 232, 463–475. [DOI](https://doi.org/10.1111/j.1365-2818.2008.02141.x). O FIRE inspira a continuação geométrica dos ramos. O código deste candidato combina esse princípio com realce Sato e análise de grafos em Python. [Créditos](NOTICE.txt).
