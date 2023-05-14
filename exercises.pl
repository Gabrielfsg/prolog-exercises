verificaSeTaNaLista(X,[X|_]).
verificaSeTaNaLista(X,[_|T]):- verificaSeTaNaLista(X,T).

maiorElemento([X|T]):- Maior is X, maiorDaLista(T,Maior).
maiorDaLista([],Maior):- write(Maior).
maiorDaLista([X|T], Maior):- X > Maior,
(
NMaior is X, maiorDaLista(T,NMaior)
) ; maiorDaLista(T,Maior).

menorElemento([X|T]):- Menor is X, menorDaLista(T,Menor).
menorDaLista([],Menor):- write(Menor).
menorDaLista([X|T],Menor):- X < Menor, NMenor is X, menorDaLista(T, NMenor) ; menorDaLista(T, Menor).

pares([],[]).
pares([X|T], Pares):- pares(T,Resto),
(
X mod 2 =:= 0, Pares = [X|Resto] ;
Pares = Resto
).

impares([],[]).
impares([X|T], Impares):- impares(T,Resto),
(
X mod 2 > 0, Impares = [X|Resto] ;
Impares = Resto
).

prefixo([],_).
prefixo([X1|T1],[X2|T2]):- X1 =:= X2, prefixo(T1,T2).

sufixo([],[]).
sufixo([X1|T1],[X2|T2]):- X1 =:= X2, sufixo(T1,T2) ; sufixo([X1|T1],T2).

todos_as([]).
todos_as([X|T]):- X = a, todos_as(T).

contem_1([X|T]):- X = 1, ! ; contem_1(T).

ultimo(X,[X]).
ultimo(X,[_|T]):- ultimo(X,T).

primeiro(X,[X|_]).

adicionaInicioLista(X,ListaEntrada,ListaSaida):- ListaSaida = [X|ListaEntrada].

adicionaFinalLista(X,_,[X]).
adicionaFinalLista(X,[Z|T],ListaSaida):- adicionaFinalLista(X,T,LS), ListaSaida = [Z|LS].

inverte([],[]).
inverte([X|T],Inv):- inverte(T,Ret), length(Ret,Tam), (Tam =:= 0, Inv = [X|Ret]  ;  append(Ret,[X], Inv)), !.

maiorFinal([X],X).
maiorFinal([X|T],Z):- maiorFinal(T,Y), (Y > X, Z  is Y ; Y < X, Z  is X).

maiorFinal2([],0).
maiorFinal2([X|T],Z):- maiorFinal2(T,Y), (X > Y, Z is X ; Y > X, Z is Y), !.

remove_duplicates([],[]).
remove_duplicates([X|T], SD):- remove_duplicates(T,R), (not(member(X,R)), SD = [X|R] ; SD = R), !.

somaDaLista([],0).
somaDaLista([X|T], Soma):- somaDaLista(T,S), Soma is S + X, !.

deletar([],Z,[]).
deletar([X|T],Z,Lista):-  deletar(T,Z,Resto), (X =:= Z, Lista = Resto ;  X \= Z, Lista = [X|Resto]).

matriz([
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0]
        ]).
        
atualizaMatriz(L,C,V,NMATRIZ):- matriz(M), upMatriz(L,C,V,M,NMATRIZ).

upMatriz(L,C,V,[X|T],NMATRIZ):-
   L =:= 1,
   atualizaColuna(C,X,V,NLinha),
   NMATRIZ = [NLinha|T].
upMatriz(L,C,V,[X|T],NMATRIZ):-
   NL is L - 1,
   upMatriz(NL,C,V,T,NM),
   NMATRIZ = [X|NM].
   
atualizaColuna(C,[X|T],V,NLinha):-
   C =:= 1,
   NLinha = [V|T].
atualizaColuna(C,[X|T],V,NLinha):-
   NC is C - 1,
   atualizaColuna(NC,T,V,N),
   NLinha = [X|N].
   
removeElementoPosicao([X|T],P,Lista):-
   P =:= 1,
   Lista = T.
removeElementoPosicao([X|T],P,Lista):-
   NP is P - 1,
   removeElementoPosicao(T,NP,NLista),
   Lista = [X|NLista].
   
   
matrizTreino([
        [5,0,1,0],
        [0,5,2,8],
        [0,2,0,3]
        ]).
        

verificaSeElementoEstaNaMatriz(E):- matrizTreino(M), verifaMatrizContem(E,M).

verifaMatrizContem(E, [X|T]):-
   verificaListaMat(E,X).
verifaMatrizContem(E, [X|T]):-
   verifaMatrizContem(E,T).

verificaListaMat(E,[E|T]).
verificaListaMat(E, [X|T]):- verificaListaMat(E,T).


somarLinha(L,S):- matrizTreino(M), somarLinhaMat(L,M,S).

somarLinhaMat(L,[X|T],S):-
   L =:= 1,
   fazSoma(X,S).
somarLinhaMat(L,[X|T],S):-
   NL is L - 1,
   somarLinhaMat(NL, T, S).

fazSoma([],0).
fazSoma([X|T],S):- fazSoma(T,Soma), S is Soma + X.


somarDuasLinhas(L1,L2, Soma):-  matrizTreino(M), somarLinha(L1,SOMA1), somarLinha(L2,SOMA2),  Soma is SOMA1 + SOMA2.


invertLinhaMatriz(L,Mat):- matrizTreino(M), invertLinha(L,M,Mat), !.

invertLinha(L,[X|T], Mat):-
   L =:= 1,
   inventerFinal(X,Nlinha),
   Mat = [Nlinha|T].
invertLinha(L,[X|T], Mat):-
   NL is L - 1,
   invertLinha(L,T, NMat),
   Mat = [X|NMat].
   
inventerFinal([],[]).
inventerFinal([X|T], Linha):- inventerFinal(T, NL), append(NL,[X],Linha).

funcaoAppend([],Lista2,ListaFinal):-  ListaFinal = Lista2.
funcaoAppend([X1|T1],Lista2,ListaFinal):- funcaoAppend(T1,Lista2,LF), ListaFinal = [X1|LF].

removePrimeiraOcorrencia(E,[E|T], T):- !.
removePrimeiraOcorrencia(E,[X|T],R):- removePrimeiraOcorrencia(E, T, NR), R = [X|NR].

comandoSelect(E,[E|T], T).
comandoSelect(E,[X|T],R):-
   comandoSelect(E, T, NR),
   R = [X|NR].
   
perm([], []).
perm(List, [X|Perm]) :-
    comandoSelect(X, List, Rest),
    perm(Rest, Perm).
    
    
multiplicarMatrizes(Mult):- matriz1(M1), matriz2(M2), fazerMult(M1,M2,Mult).

fazerMult([],[],[]).
fazerMult([X1|T1],[X2|T2], M):-
   fazerMult(T1,T2,NM),
   multDuasLinhas(X1,X2,NL),
   M = [NL|NM].

multDuasLinhas([],[],[]).
multDuasLinhas([X1|T1],[X2|T2],L):- multDuasLinhas(T1,T2,NL), NX is X1 * X2, L = [NX|NL].


matriz1([
        [1,4,7],
        [2,5,8],
        [3,6,9]
        ]).

matriz2([
        [1,2,3],
        [4,5,6],
        [7,8,9]
        ]).

somarMatriz(Soma):- matriz1(M1), matriz2(M2), fazerSoma(M1,M2,Soma).

fazerSoma([],[],[]).
fazerSoma([L1|T1],[L2|T2],Soma):-
   fazerSoma(T1,T2, NS),
   somaLinhas(L1,L2,R),
   Soma = [R|NS].

somaLinhas([],[],[]).
somaLinhas([X1|T1],[X2|T2],R):- somaLinhas(T1,T2,NR), X3 is X1 + X2, R = [X3|NR].

constroiMatriz(0,C,[]).
constroiMatriz(L,C,M):-
   LT is L - 1,
   constroiMatriz(LT,C,NM),
   constroiColunas(C, NC),
   M = [NC|NM], !.

constroiColunas(0,[]).
constroiColunas(C,L):- CT is C - 1, constroiColunas(CT,NL), L = [0|NL].


upMatNew(L,C,V,Resp):- matriz2(M), updateMat(L,C,V,M,Resp).
   
updateMat(L,C,V,[X|T],Resp):-
   L =:= 1,
   atualizarAMat(C,V, X, NL),
   Resp = [NL|T].
updateMat(L,C,V,[X|T],Resp):-
   NL is L - 1,
   updateMat(NL,C,V,T,NResp),
   Resp = [X|NResp], !.

atualizarAMat(0,V,T,T).
atualizarAMat(C,V,[X|T],Ret):- CT is C - 1, atualizarAMat(CT,V,T,NRET), (CT =:= 0 , Ret = [V|NRET] ; Ret = [X|NRET]).


mat1([
       [1,2,3],
       [4,8,6],
       [27,8,9]
       ]).
       
maiorPorLinha(L):- mat1(M), buscaMaiorPorLinha(M,L).

buscaMaiorPorLinha([],[]).
buscaMaiorPorLinha([X|T],L):-
   buscaMaiorPorLinha(T,NL),
   maiorLinha(X,E),
   L = [E|NL].

  
maiorLinha([],0).
maiorLinha([X|T],E):- maiorLinha(T,NE), (X > NE, E is X ; E is NE), !.

diagonalMat(R):- mat1(M), length(M,Var), retornaDiagonal(M,0,R).

retornaDiagonal([],V,[]).
retornaDiagonal([X|T],V,R):-
   NV is V + 1,
   retornaDiagonal(T,NV,NR),
   retornaElementoDia(X,NV,E),
   R = [E|NR].


retornaElementoDia([X|T], 1, X).
retornaElementoDia([X|T], V, E):- NV is V - 1, retornaElementoDia(T, NV, E).


grafo(1,2,1.5).
grafo(1,3,2.0).
grafo(2,3,1.1).
grafo(3,4,3.2).
grafo(5,6,4.1).
grafo(4,5,3.2).
grafo(3,5,3.2).
grafo(6,7,2.1).
grafo(7,8,3.7).
grafo(8,6,1.0).
grafo(8,9,2.8).

vizinhos(V1,V2,D):- grafo(V1,V2,D) ; grafo(V2,V1,D).

conectados(V1,V2,Visi):-
   vizinhos(V1,V2,_).
conectados(V1,V2,Visi):-
   grafo(V1,Y,_),
   not(member(Y,Visi)),
   conectados(Y,V2,[Y|Visi]).
   
   
percurso(V1,V2,Visitados,Perc,Dist):-
   vizinhos(V1,V2,D),
   append(Perc,[V1,V2],NPERC),
   NDIST is Dist + D,
   write(NPERC), nl,
   write(NDIST).
percurso(V1,V2, Visitados, Perc,Dist):-
   grafo(V1,Y,D),
   not(member(Y,Visitados)),
   append(Perc,[V1],NPERC),
   NDIST is Dist + D,
   percurso(Y,V2, [Y|Visi], NPERC,NDIST).
   
   
animal(urso).
animal(guaxinim).
animal(peixe).
animal(peixinho).
animal(lince).
animal(raposa).
animal(coelho).
animal(veado).

come(peixe,peixinho).
come(peixinho,alga).
come(guaxinim,peixe).
come(urso,peixe).
come(raposa,coelho).
come(lince,veado).
come(coelho,grama).
come(urso,raposa).
come(veado,grama).
come(urso,veado).
come(urso,guaxinim).

planta(alga).
planta(grama).

presa(X):- animal(Y), come(Y,X).


aluno(joao, calculo).
aluno(maria, calculo).
aluno(joel, programacao).
aluno(joel, estrutura).

frequenta(joao, puc).
frequenta(maria, puc).
frequenta(joel, ufrj).

professor(carlos, calculo).
professor(ana_paula, estrutura).
professor(pedro, programacao).

funcionario(pedro, ufrj).
funcionario(ana_paula, puc).
funcionario(carlos, puc).


alunosProfessor(Prof, Aluno):- aluno(Aluno, Materia), professor(Prof, Materia), frequenta(Aluno,  Facu), funcionario(Prof, Facu).

potencia(N,2,R):- R is N * N, !.
potencia(N,P,R):- NP is P - 1, potencia(N,NP,NR),  R is NR * N, !.

montaIntervalo(X,NX,[]):- NX is X - 1.
montaIntervalo(X,Y,I):- NY is Y - 1, montaIntervalo(X,NY,NI), append(NI,[Y],I), !.

elementoPerIndice(0,[X|T],X).
elementoPerIndice(I,[X|T],E):- NI is I - 1, elementoPerIndice(NI,T,E).

progenitor(maria, jose).
progenitor(joao, jose).
progenitor(joao, ana).
progenitor(jose, julia).
progenitor(jose, iris).
progenitor(iris, jorge).

filho(X,Y):- progenitor(X,Y).



conectado(1,2).
conectado(3,4).
conectado(5,6).
conectado(7,8).
conectado(9,10).
conectado(12,13).
conectado(13,14).
conectado(15,16).
conectado(17,18).
conectado(19,20).
conectado(4,1).
conectado(6,3).
conectado(4,7).
conectado(6,11).
conectado(14,9).
conectado(11,15).
conectado(16,12).
conectado(14,17).
conectado(16,19).

caminho(X,Y):-
   conectado(X,Y).
caminho(X,Y):-
   conectado(X,Z),
   caminho(Z,Y).

removeMultiplos([],D,[]).
removeMultiplos([X|T],D,L):- removeMultiplos(T,D,NL), (Resto is X mod D, Resto \= 0, L = [X|NL]  ; L = NL), !.

removePrimos([],[]).
removePrimos([X|T],L):- removePrimos(T,NL), (primo(X), L = NL ; L = [X|NL]), !.

intervaloPrimosUm(V,Lista):- geraSeq(V,Seq), intervaloPrimosDois(Seq,Lista).

intervaloPrimosDois([],[]).
intervaloPrimosDois([X|T], Lista):- intervaloPrimosDois(T,NLista), (primo(X), Lista = [X|NLista] ; Lista = NLista), !.

primo(1).
primo(N):- N > 1, NN is N - 1, geraSeq(NN,L), verPrimo(N,L).

verPrimo(N,[]).
verPrimo(N,[X|T]):- Resto is N mod X, Resto \= 0, verPrimo(N,T).

geraSeq(1,[]).
geraSeq(N,L):- Nn is N - 1, geraSeq(Nn,NL), L = [N|NL], !.

homem(pedro).
homem(marcos).
homem(ze).
mulher(maria).
mulher(joana).
idade(ze,30).
idade(maria,40).
idade(marcos,20).
idade(pedro,25).
idade(joana,28).
gosta(ze,aventura).
gosta(maria,comedia).
gosta(joana,romance).
gosta(marcos,terror).
gosta(marcos,romance).
gosta(pedro,romance).
gosta(maria,romance).
opcao(ze,20,40).
opcao(pedro,25,55).
opcao(marcos,20,21).
opcao(maria,25,55).
opcao(joana,28,30).

afinidadeFilme(X,Y):- gosta(X,Z), gosta(Y,Z).
casal(X,Y):- homem(X), mulher(Y) ; homem(Y), mulher(X).

elementosRepetidos(Lista):- temRepetidos(Lista,[]).

temRepetidos([],Lista).
temRepetidos([X|T],Lista):- nAopertence(X,Lista), NLista =  [X|Lista], temRepetidos(T,NLista).

nAopertence(E,[]).
nAopertence(E,[X|T]):- E \= X, nAopertence(E,T).

tamanho([],0).
tamanho([X|T],Tam):-  tamanho(T,NTam), Tam is  NTam + 1.

par(Num):- Resto is Num mod 2, Resto =:= 0.

meioLista(Lista,E):- tamanho(Lista,Tam), (not(par(Tam)), I is (Tam / 2) - 0.5, buscaMeioL(Lista,I,E) ; par(Tam), I is (Tam / 2), buscaMeioL(Lista,I,E)).

buscaMeioL([X|T],0.0,X).
buscaMeioL([X|T],0,X).
buscaMeioL([X|T],I,E):- NI is I - 1, buscaMeioL(T,NI,E).

deleta1(X,[X|T],T).
deleta1(E,[X|T],Lista):- deleta1(E,T,NLista), Lista = [X|NLista], !.

deletaAllOc(E,[],[]).
deletaAllOc(E,[X|T],Lista):- deletaAllOc(E,T,Nlista), (E =:= X, Lista = Nlista ; Lista = [X|Nlista]), !.

subsPT(E,X,[X|T],[E|T]).
subsPT(E,EL,[X|T],R):- subsPT(E,EL,T,NR), R = [X|NR].

mediaListaForne(Lista,M):- somaDaListaTotal(Lista,Soma),  tamanho(Lista,Tam), M is Soma/Tam.

somaDaListaTotal([], 0).
somaDaListaTotal([X|T],Soma):- somaDaListaTotal(T, Nsoma), Soma is Nsoma + X.

mediana(Lista,Mediana):- ordemCrescente(Lista,Crescente), calcMediana(Crescente,Mediana), !.

calcMediana(Crescente,Mediana):- length(Crescente, Tam), Resto is Tam mod 2, (Resto =:= 0, I1 is (Tam / 2) - 1, I2 is (Tam / 2), medianaPar(Crescente,I1,I2,Mediana) ; I is (Tam / 2) - 0.5, buscaMeioL(Crescente,I,Mediana)).

medianaPar(L,I1,I2,M):- indiceRetEle(L,I1,E1), indiceRetEle(L,I2,E2), M is (E1+E2)/2.

indiceRetEle([X|T],0, X).
indiceRetEle([X|T],I,E):- NI is I - 1, indiceRetEle(T,NI,E).

ordemCrescente([],[]).
ordemCrescente([X|T],Lista):- ordemCrescente(T,NL), moveListaCres(X,NL,RNL),Lista = RNL, !.

moveListaCres(E,[],[E]).
moveListaCres(E,[X|T],RNL):- ( E < X, append([E],[X|T],RNL) ; moveListaCres(E,T,NRNL), RNL  = [X|NRNL]).

% Problemas comuns vivenciados em uma assistencia técnica
problema('computador_nao_liga').
problema('computador_apitando').
problema('monitor_nao_da_video').
problema('teclado_nao_conecta').
problema('teclado_funcionando_mal').
problema('mouse_nao_conecta').
problema('mouse_nao_clica').
problema('impressora_nao_imprime').
problema('hora_do_computador_desregulada').
problema('computador_sem_internet').
problema('super_aquecimento_do_computador').
problema('lentidao_no_computador').
problema('tela_azul_da_morte').
problema('falha_no_disco_rigido_ou_SSD').
problema('falta_de_espaco_no_disco_rigido').
problema('virus_ou_malware_no_computador').
problema('som_nao_funcionando_no_computador').
problema('defeito_monitor').
problema('reinicializacao_inesperada_do_sistema').

/* ----------------------------------------------------------------------- */

motivoProblema('computador_nao_liga',Prob,Sol) :- problema_bios(Prob), solucao_problema(Prob,Sol).
motivoProblema('computador_nao_liga',Prob,Sol) :- problema_ram(Prob), solucao_problema(Prob,Sol).

motivoProblema('computador_nao_liga',Prob,Sol) :- problema_bateria_mobo(Prob), solucao_problema(Prob,Sol).
motivoProblema('computador_nao_liga',Prob,Sol) :- problema_processador(Prob), solucao_problema(Prob,Sol).
motivoProblema('computador_nao_liga',Prob,Sol) :- problema_fonte('computador nao liga', Prob), solucao_problema(Prob,Sol).

motivoProblema('computador_apitando',Prob,Sol) :- problema_ram(Prob), solucao_problema(Prob,Sol).

motivoProblema('computador_apitando',Prob,Sol) :- problema_bios(Prob), solucao_problema(Prob,Sol).
motivoProblema('computador_apitando',Prob,Sol) :- problema_bateria_mobo(Prob), solucao_problema(Prob,Sol).
motivoProblema('computador_apitando',Prob,Sol) :- problema_processador(Prob), solucao_problema(Prob,Sol).

motivoProblema('monitor_nao_da_video',Prob,Sol) :- problema_ram(Prob),solucao_problema(Prob,Sol).
motivoProblema('monitor_nao_da_video',Prob,Sol) :- problema_gpu(Prob),solucao_problema(Prob,Sol).
motivoProblema('monitor_nao_da_video',Prob,Sol) :- problema_monitor(Prob), solucao_problema(Prob,Sol).

motivoProblema('teclado_nao_conecta',Prob,Sol) :- problema_periferico(Prob), solucao_problema(Prob,Sol).
motivoProblema('teclado_nao_conecta',Prob,Sol) :- problema_teclado(Prob), solucao_problema(Prob,Sol).
motivoProblema('teclado_funcionando_mal',Prob,Sol) :- problema_periferico(Prob), solucao_problema(Prob,Sol).
motivoProblema('teclado_funcionando_mal',Prob,Sol) :- problema_teclado(Prob), solucao_problema(Prob,Sol).

motivoProblema('mouse_nao_conecta',Prob,Sol) :- problema_periferico(Prob), solucao_problema(Prob,Sol).
motivoProblema('mouse_nao_conecta',Prob,Sol) :- problema_mouse(Prob), solucao_problema(Prob,Sol).
motivoProblema('mouse_nao_clica',Prob,Sol) :- problema_periferico(Prob), solucao_problema(Prob,Sol).
motivoProblema('mouse_nao_clica',Prob,Sol) :- problema_mouse(Prob), solucao_problema(Prob,Sol).

motivoProblema('impressora_nao_imprime',Prob,Sol) :- problema_impressora(Prob), solucao_problema(Prob,Sol).

motivoProblema('hora_do_computador_desregulada',Prob,Sol) :- problema_computador_hora(Prob), solucao_problema(Prob,Sol).

motivoProblema('computador_sem_internet',Prob,Sol) :- problema_computador_internet(Prob), solucao_problema(Prob,Sol).

motivoProblema('super_aquecimento_do_computador',Prob,Sol) :- problema_fonteInsuficiente(Prob), solucao_problema(Prob,Sol).

motivoProblema('lentidao_no_computador',Prob,Sol) :- problema_lentidao(Prob), solucao_problema(Prob,Sol).

motivoProblema('tela_azul_da_morte',Prob,Sol) :- problema_telaAzul(Prob), solucao_problema(Prob,Sol).

motivoProblema('falha_no_disco_rigido_ou_SSD',Prob,Sol) :- problema_hdssd(Prob), solucao_problema(Prob,Sol).

motivoProblema('falta_de_espaco_no_disco_rigido',Prob,Sol) :- problema_hdssd_cheio(Prob), solucao_problema(Prob,Sol).

motivoProblema('virus_ou_malware_no_computador',Prob,Sol) :- problema_hdssd(Prob), solucao_problema(Prob,Sol).

motivoProblema('som_nao_funcionando_no_computador',Prob,Sol) :- problema_som(Prob), solucao_problema(Prob,Sol).

motivoProblema('defeito_monitor',Prob,Sol) :- problema_tela(Prob), solucao_problema(Prob,Sol).

motivoProblema('reinicializacao_inesperada_do_sistema',Prob,Sol) :- problema_ram(Prob), solucao_problema(Prob,Sol).
motivoProblema('reinicializacao_inesperada_do_sistema',Prob,Sol) :- problema_processador(Prob), solucao_problema(Prob,Sol).
motivoProblema('reinicializacao_inesperada_do_sistema',Prob,Sol) :- problema_fonte('computador nao liga', Prob), solucao_problema(Prob,Sol).

/* ----------------------------------------------------------------------- */
% Fatos que descrevem os possiveis problemas para cada componente da maquina
problema_ram('memoria_ram_com_defeito').
problema_ram('memoria_ram_desencaixada').
problema_ram('memoria_ram_suja').
problema_ram('slot_memoria_ram_com_defeito').

problema_bios('bios_pode_estar_desconfigurada').
problema_bios('bios_pode_estar_corrompida').

problema_bateria_mobo('bateria_descarregada').

problema_processador('processador_com_defeito').
problema_processador('processador_desencaixado_do_socket').
problema_processador('ventoinha_desencaixada_do_processador').

problema_fonte('fonte_de_alimentação_em_curto').
problema_fonte('computador_nao_liga', 'fonte_nao_fornece_energia_suficiente_aos_componentes').

problema_mobo('problema_no_socket_do_processador').
problema_mobo('pino_do_socket_torto').
problema_mobo('problema_no(s)_chipset(s)').

problema_gpu('gpu_nao_esta_encaixada_corretamente').
problema_gpu('gpu_esta_suja').
problema_gpu('gpu_com_defeito_eletronico').

problema_monitor('cabo_de_video_com_defeito').
problema_monitor('monitor_com_defeito_eletronico').

problema_periferico('dispositivo_em_curto').
problema_periferico('usb_ou_p2_da_placa-mae_com_problema').
problema_periferico('usb_ou_p2_do_dispositivo_com_problema').
problema_periferico('driver_com_problemas').
problema_periferico('dispositivo_com_defeito_eletronico').
problema_periferico('algum_outro_periferico_na_maquina_em_curto_esta_impedindo_o_funcionamento_correto').

problema_teclado('teclas_oxidadas_ou_sujas').
problema_teclado('teclado_mecanico_switch_com_mau_contato_ou_queimado').
problema_teclado('teclado_membrana_trilha_da_membrana_pode_estar_rompida').

problema_mouse('botoes_oxidados_ou_sujos').
problema_mouse('botoes_queimados').

problema_impressora('papel_agarrado_ou_obstruído').
problema_impressora('cartucho_entupido_ou_vazio').
problema_impressora('cartucho_queimado').
problema_impressora('spooler_de_impressão_com_problemas').

problema_computador_hora('bateria_da_placa_mae_descarregada').

problema_computador_internet('cabo_com_problema_ou_sem_sinal_wi_fi').
problema_computador_internet('provedor_sem_internet').
problema_computador_internet(X) :- problema_roteador(X).

problema_roteador('roteador_com_as_configurações_de_rede_incorretas').
problema_roteador('roteador_com_defeito_na antena').
problema_roteador('roteador_com_defeito_eletronico').

problema_fonteInsuficiente('fonte_com_poder_insuficiente').

problema_lentidao('computador_lento').

problema_telaAzul('computador_ou_notebook_com_tela_azul_da_morte').

problema_hdssd('hd_ou_sdd_lento').
problema_hdssd('hd_ou_sdd_nao_foram_identificados_pelo_SO').

problema_hdssd_cheio('hd_ou_sdd_lotado').

problema_virus('computador_com_virus').

problema_som('caixas_de_som_nao_funciona').
problema_som('fones_de_ouvido_nao_funciona').

problema_tela('tela_piscando_ou_tremendo').

/* -------------------------------------------------------------------- */
% Soluções para os problemas

solucao_problema('memoria_ram_com_defeito',Sol) :- Sol ='troque a memoria ram'.
solucao_problema('memoria_ram_com_defeito',Sol) :- Sol ='limpe a memoria ram'.
solucao_problema('memoria_ram_desencaixada',Sol) :- Sol ='retira e encaixe novamente a memoria'.
solucao_problema('memoria_ram_suja',Sol) :- Sol ='limpe as memorias com limpa contato'.
solucao_problema('slot_memoria_ram_com_defeito',Sol) :- Sol ='Troque o slot da placa-mae ou troque de placa-mae'.
solucao_problema('bios_pode_estar_desconfigurada',Sol) :- Sol ='Troque a bateria da placa mae e reconfigure a BIOS'.
solucao_problema('bios_pode_estar_corrompida',Sol) :- Sol ='Resete a BIOS'.
solucao_problema('bateria_descarregada',Sol) :- Sol ='Troque a bateria'.
solucao_problema('processador_com_defeito',Sol) :- Sol ='Troque o processador'.
solucao_problema('processador_desencaixado_do_socket',Sol) :- Sol ='Retire e encaixe o processador no socket corretamente.'.
solucao_problema('ventoinha_desencaixada_do_processador',Sol) :- Sol ='Retire e encaixe a ventoinha corretamente'.
solucao_problema('fonte_de_alimentação_em_curto',Sol) :- Sol ='Troque a fonte de alimentacao'.
solucao_problema('problema_no_socket_do_processador',Sol) :- Sol ='Realize a troca do socket da placa-mae ou troque de placa-mae'.
solucao_problema('pino_do_socket_torto',Sol) :- Sol ='Desentorte os pinos cuidadosamente'.
solucao_problema('problema_no(s)_chipset(s)',Sol) :- Sol ='Troque de placa-mae'.
solucao_problema('gpu_nao_esta_encaixada_corretamente',Sol) :- Sol ='Retire e encaixe novamente a gpu'.
solucao_problema('gpu_esta_suja',Sol) :- Sol ='Limpe a GPU com limpa-contato e/ou alcool isopropilico'.
solucao_problema('gpu_com_defeito_eletronico',Sol) :- Sol ='Troque de GPU'.
solucao_problema('cabo_de_video_com_defeito',Sol) :- Sol ='Troque de cabo'.
solucao_problema('monitor_com_defeito_eletronico',Sol) :- Sol ='Troque de monitor'.
solucao_problema('dispositivo_em_curto',Sol) :- Sol ='Troque o dispositivo'.
solucao_problema('algum_periferico_na_maquina_em_curto_esta_impedindo_o_funcionamento_correto',Sol) :- Sol ='Verifique cada periferico, testando-os individualmente na maquina'.
solucao_problema('driver_com_problemas',Sol) :- Sol ='Reinstale o driver correto'.
solucao_problema('usb_ou_p2_da_placa-mae_com_problema',Sol) :- Sol ='Teste o dispositivo em outras portas'.
solucao_problema('teclas_oxidadas_ou_sujas',Sol) :- Sol ='Desmonte o teclado e limpe as teclas com alcool isopropilico'.
solucao_problema('teclado_mecanico_switch_com_mau_contato_ou_queimado',Sol) :- Sol ='Limpe os switches com limpa contato'.
solucao_problema('teclado_membrana_trilha_da_membrana_pode_estar_rompida',Sol) :- Sol ='Faca a reconstrucao da trilha ou troque a membrana'.
solucao_problema('botoes_oxidados_ou_sujos',Sol) :- Sol ='Desmonte o mouse e limpe os botoes com alcool isopropilico'.
solucao_problema('botoes_queimados',Sol) :- Sol ='Troque o mouse'.
solucao_problema('dispositivo_com_defeito_eletronico',Sol) :- Sol ='Troque o dispositivo'.
solucao_problema('papel_agarrado_ou_obstruído',Sol) :- Sol ='Abra a tampa da impressora e retire o papel'.
solucao_problema('cartucho_entupido_ou_vazio',Sol) :- Sol ='Recarregue ou desentupa o cartucho'.
solucao_problema('cartucho_queimado',Sol) :- Sol ='Troque o cartucho'.
solucao_problema('spooler_de_impressão_com_problemas',Sol) :- Sol ='Reinicie o spooler de impressao do sistema operacional'.
solucao_problema('spooler_de_impressão_com_problemas',Sol) :- Sol ='Reinicie o spooler de impressao do sistema operacional'.
solucao_problema('bateria_da_placa-mae_descarregada',Sol) :- Sol ='Troque a bateria'.
solucao_problema('cabo_com_problema_ou_sem_sinal_wi_fi',Sol) :- Sol ='Verifique a integridade do cabo ou sinal do wi-fi com outro dispositivo'.
solucao_problema('provedor_sem_internet',Sol) :- Sol ='Contate o seu provedor'.
solucao_problema('roteador_com_as_configurações_de_rede_incorretas',Sol) :- Sol ='Realize as configuracoes do roteador baseado nas utilizadas pelo provedor'.
solucao_problema('roteador_com_defeito_na_antena',Sol) :- Sol ='Troque a antena ou o roteador'.
solucao_problema('roteador_com_defeito_eletronico',Sol) :- Sol = 'Troque o roteador'.
solucao_problema('fonte_nao_fornece_energia_suficiente_aos_componentes',Sol) :- Sol = 'Troque a fonte de alimentacao'.
solucao_problema('fonte_com_poder_insuficiente',Sol) :- Sol = 'Compre outra fonte'.
solucao_problema('computador_lento',Sol) :- Sol = 'Finalize os processos nao utilizados'.
solucao_problema('computador_ou_notebook_com_tela_azul_da_morte',Sol) :- Sol = 'Reinicie o computador'.
solucao_problema('hd_ou_sdd_lento',Sol) :- Sol = 'Desfragmente seu HD OU SSD'.
solucao_problema('hd_ou_sdd_nao_foram_identificados_pelo_SO',Sol) :- Sol = 'Verifique se os modulos de armazenamento foram encaixados corretamente'.
solucao_problema('hd_ou_sdd_lotado',Sol) :- Sol = 'Remova dados nao utilizados no HD OU SSD'.
solucao_problema('computador_com_virus',Sol) :- Sol = 'Instale um antivirus ou formate o computador'.
solucao_problema('caixas_de_som_nao_funciona',Sol) :- Sol = 'Verifique se as caixas de som estao conectados corretamente, ou se os drivers de som estão instalados'.
solucao_problema('fones_de_ouvido_nao_funciona',Sol) :- Sol = 'Verifique se o fone de ouvido esta conectado corretamente, ou se os drivers do fone/som estão instalados'.
solucao_problema('tela_piscando_ou_tremendo',Sol) :- Sol = 'Troque o monitor, ou verifique se o monitor esta com os leds estao funcionando, ou se ele tem algum outro componente estragado.'.


/* -------------------------------------------------------------------- */

consultarProblema :-
                  write(''),
                  write('Digite o problema que o seu computador possui.'), nl,
                  write('Exemplo: computador nao liga. '), nl,
                  write('O problema deve ter "." no final'), nl,
                  write('Para finalizar a entrada de dados, entre com: "f."'), nl,
                  write('Para sair do programa: "sair."'), nl,
                  write('Problemas: '), nl,
                  listarProblemas,
                  lerProblemas([]).

listarProblemas :-
    problema(X),
    write(X), nl,
    fail.
listarProblemas.

lerProblemas(Problemas):-
read(Problema), identificarProblemas(Problema,Problemas).

identificarProblemas(sair,_).

identificarProblemas(f,Problemas):-
  (validaTamanhoLista(Problemas),buscarProblemas(Problemas) ; write('Deve ter ao menos um problema na lista. '), nl, lerProblemas(Problemas)).

identificarProblemas(Problema,Problemas):-
  problema(Problema), (not(pertence(Problema,Problemas)), write('Problema ja foi adicionado. '), nl, lerProblemas(Problemas) ; lerProblemas([Problema|Problemas])).

identificarProblemas(Problema,Problemas):-
  validaFatoComParametro(problema(Problema)),
  write('O problema desconhecido na base de dados. '), nl,
  lerProblemas(Problemas).

validaFatoComParametro(Predicado) :-
    Predicado, !, fail.
validaFatoComParametro(_).

validaTamanhoLista(Lista):-
  length(Lista,Tam), Tam > 0.

buscarProblemas(Entrada):-
  preencheLista(Entrada,Lista),
  write('Problemas:'), nl,
  listaProblemasSecundarios(Lista),
  write('Digite os problemas de acordo com os listados com "." no final.'), nl,
  write('Digite "f." para finalizar a entrada.'), nl,
  write('Digite "sair." para finalizar o programa.'), nl,
  lerProblemasSecundarios([],Lista).

preencheLista([], []).
preencheLista([X|T], ListaResultado) :-
    preencheLista(T, ListaTemp),
    findall(Prob, motivoProblema(X, Prob, _), ListaProblemas),
    append(ListaProblemas, ListaTemp, ListaResultado).

listaProblemasSecundarios([X|T]):- length([X|T],Tamanho),
                                    Tamanho > 0, write(X),
                                    nl, listaProblemasSecundarios(T) ;  !.

lerProblemasSecundarios(ListaEntrada,ListaProblemas):-
   read(Problema), verificaSeEProblema(Problema, ListaEntrada,ListaProblemas).

verificaSeEProblema(sair,_,_).
verificaSeEProblema(f,ListaEntrada,ListaProblemas):-
   (validaTamanhoLista(ListaEntrada), buscarSolucoes(ListaEntrada) ; write('Deve ter ao menos um problema na lista. '), nl, lerProblemasSecundarios(ListaEntrada,ListaProblemas)).

verificaSeEProblema(Problema,ListaEntrada,ListaProblemas):-
   member(Problema, ListaProblemas),
   (not(pertence(Problema,ListaEntrada)), write('Problema ja foi adicionado. '), nl, lerProblemasSecundarios(ListaEntrada,ListaProblemas) ; lerProblemasSecundarios([Problema|ListaEntrada],ListaProblemas)).

verificaSeEProblema(Problema,ListaEntrada,ListaProblemas):-
   not(member(Problema, ListaProblemas)),
   write('O problema desconhecido na base de dados. '), nl,
   lerProblemasSecundarios(ListaEntrada,ListaProblemas).

buscarSolucoes(ListaProblemas) :-
    member(Problema, ListaProblemas),
    solucao_problema(Problema, Solucao),
    write('Problema: '), write(Problema), write(' Solucao: '), write(Solucao), nl,
    fail.
buscarSolucoes(_):- consultarProblema.

pertence(_, []) :- !.
pertence(X, [X|_]) :- !, fail.
pertence(X, [_|T]) :- pertence(X, T).


criaTabuleiro(N,Tab) :- N > 0, criaTabuleiro(N,N,[],Tab).
criaTabuleiro(0,_,Acc,Tab) :- Tab = Acc, !.
criaTabuleiro(X,N,Acc,Tab) :- length(L,N), maplist(=(0),L), X1 is X - 1, criaTabuleiro(X1,N,[L|Acc],Tab).

imprimeTab(Tab) :-
   imprimeLinha(Tab).

imprimeLinha([]).
imprimeLinha([L|[]]) :-
   write(L).

imprimeLinha([L|R]) :-
   write(L),
   nl,
   imprimeLinha(R).

programaTeste :- criaTabuleiro(5,Tab),
            saltosCavalo(Tab,2,2).

programa :- write('Defina o tamanho do tabuleiro: '),
            read(N),
            criaTabuleiro(N,Tab), 
            write('Tabela criada: '),
            nl,
            imprimeTab(Tab),
            nl,
            write('Defina a linha inicial (começa em 0): '),
            read(Lin),
            write('Defina a coluna inicial (começa em 0): '),
            read(Col),
            nl,
            write('Saltando...'),
            saltosCavalo(Tab,Lin,Col).

movimento((-1,-2)).
movimento((1,-2)).
movimento((-2,-1)).
movimento((-2,1)).
movimento((-1,2)).
movimento((1,2)).
movimento((2,-1)).
movimento((2,1)).

saltosCavalo(Tab,Lin,Col) :- length(Tab,N),
                             atualiza((Lin,Col),1,Tab,NewTab),
                             Visitados = [(Lin,Col)],
                             prepararSalto(N,NewTab,Lin,Col,2,Visitados).

prepararSalto(N,Tab,Lin,Col,Saltos,Visitados) :- Saltos =:= N*N + 1,!,write('Resultado apos os saltos: '),nl,imprimeTab(Tab);
                                                 faz_o_L(N,Lin,Col,Tab,Saltos,NewTab,NewLin,NewCol,Visitados,NV),
                                                 NewSaltos is Saltos + 1,
                                                 prepararSalto(N,NewTab,NewLin,NewCol,NewSaltos,NV).
                                                 
faz_o_L(Tam,Lin,Col,Tab,Salto,NewTab,NewLin,NewCol,Visitados,NV) :- movimento((M5,M6)),
                                                   NewLin is Lin + M5, NewLin >= 0, NewLin < Tam,
                                                   NewCol is Col + M6, NewCol >= 0, NewCol < Tam,
                                                   \+ member((NewLin,NewCol),Visitados), 
                                                   preencheRastro((NewLin,NewCol),Tab,Salto,NewTab,Visitados,NV).

preencheRastro((M5,M6),Tab,S,NewTab,Visitados,VisitadosRastro) :-
                                                                 atualiza((M5,M6),S,Tab,NewTab),
                                                                 append(Visitados,[(M5,M6)],VisitadosRastro).
                                                                 
atualiza((Lin,Col),N,Tab,NovaTab) :-
                                   nth0(Lin,Tab,OldLin),
                                   atualizaCol(OldLin,Col,N,NewLin),
                                   atualizaLin(Tab,Lin,NewLin,NovaTab).
atualizaCol(Lin,Col,N,NewLin) :-
                              nth0(Col,Lin,_,R),
                              nth0(Col,NewLin,N,R).

atualizaLin(Tab,Lin,NewLin,NewTab) :-
                                    nth0(Lin,Tab,_,R),
                                    nth0(Lin,NewTab,NewLin,R).
