/**
 * Autores: Geraldo José Ferreira Chagas Júnior
 *
 *
 * Data: 09/2016
 *
 * Programa de Algorítmo Genetico para resolver tsp
 * Necessário instalação da biblioteca libxml
 * Instalação da biblioteca no debian
 *           $> sudo apt-get install libxml
 *
 * compilação:
 *           $> g++ `xml2-config --cflags --libs` -o tsp tsp.cpp
 *
 * 
 * Testes:
 *      1) Mutação 2opt, apenas se melhorar a rota X Mutção 2opt sempre 
 *                                                   (Sempre multacção no melhor X aleatório para o melhor)
 *      2) Nova geração substituindo os piores     X Nova geração substituindo a menor diferença
 *      3) Sorteio de pares, todos o mesmo peso    X Sorteio de pares com peso
 *      4) cross over de pares aleatórios          X Sorteio de pares p1 e p2 em ordem
 *                                                   (melhor prmeiro X pior primeiro)
 */
#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <string.h>
#include <time.h>

#ifdef LIBXML_TREE_ENABLED

int compare(const void *x, const void *y);

/*******************************************************
classe de TMapaGenes. Todas as distâncias entre os genes
********************************************************/
class TMapaGenes
{
   public: 
      static const double infinito = 1000000.00;

   private:
      double **VP_mapaDist;
      int VP_qtdeGenes;

   //Metodos Privados
   int getNumGeneDoArquivo(xmlDocPtr doc, xmlNode * a_node)
   {
      xmlNode *cur_node = NULL;
      xmlChar *key;
      int val;

      srand(time(NULL));
      for (cur_node = a_node; cur_node; cur_node = cur_node->next) 
      {
         if (cur_node->type == XML_ELEMENT_NODE) 
         {
            if (!xmlStrcmp(cur_node->name, (xmlChar *)"description"))
            {
               key = xmlNodeListGetString(doc, cur_node->xmlChildrenNode, 1);
               val = atoi((char *)key);
               xmlFree(key); 
               return val;
            }
         }
      }
      return 0;
   }

   void preencheMapaDist (int geneOri, xmlDocPtr doc, xmlNode * a_node)
   {
      xmlNode *cur_node = NULL;
      xmlChar *uri;
      xmlChar *key;
      char dist[25];
      int geneDest;
      float df;

      for (cur_node = a_node; cur_node; cur_node = cur_node->next) 
      {
         if (cur_node->type == XML_ELEMENT_NODE)
         {
            uri = xmlGetProp(cur_node, (const xmlChar *)"cost");
            sprintf(dist, "0%s", uri);
            df = atof (dist);

            key = xmlNodeListGetString(doc, cur_node->xmlChildrenNode, 1);
            geneDest = atoi((char *)key);

            set_distancia(geneOri, geneDest, df);

            xmlFree(key); 
	    xmlFree(uri);
         }

      }
   }

   void preencheMapa(xmlDocPtr doc, xmlNode * a_node)
   {
      xmlNode *cur_node = NULL;
      int gene=0;
      int dist;

      for (cur_node = a_node; cur_node; cur_node = cur_node->next) 
      {
         if (cur_node->type == XML_ELEMENT_NODE) 
         {
            if (!xmlStrcmp(cur_node->name, (xmlChar *)"vertex"))
            {
               preencheMapaDist(gene++, doc, cur_node->children);               
            }
         }
         preencheMapa(doc, cur_node->children);
      }
   }

   public:
      TMapaGenes ()
      {
        VP_qtdeGenes = -1;
      }

      TMapaGenes (int numGenes)
      {
         inicializa (numGenes);
      }

      int get_qtdeGenes () { return VP_qtdeGenes; };

      void carregaDoArquivo(char *nomeArquivo)
      {
         int numGenes;
         xmlDoc *doc = NULL;
         xmlNode *root_element = NULL;

         // Lendo o arquivo 
         doc = xmlReadFile(nomeArquivo, NULL, 0);

         if (doc == NULL) 
         {
            printf("Erro ao carregar o arquivo %s\n", nomeArquivo);
            return;
         }

         // Obtendo o elemento root
         root_element = xmlDocGetRootElement(doc);
         
         //Obtendo o número de genes no arquivo
         VP_qtdeGenes = getNumGeneDoArquivo(doc, root_element->children);
          
         //Alocando a tabela
         inicializa (VP_qtdeGenes);

         //preenchendo a tabela com os valores da distáncia
         preencheMapa(doc, root_element->children);

         //liberando documento
         xmlFreeDoc(doc);
         // liberando as variaveis lobais
         xmlCleanupParser();

      }

      void inicializa (int numGenes) 
      {
         int i;
         int j;
         VP_qtdeGenes = numGenes;

         VP_mapaDist = (double **) malloc(numGenes*sizeof(double *)); 
         for (i=0; i<VP_qtdeGenes; i++)
         {
            VP_mapaDist [i] = (double *) malloc(numGenes*sizeof(double));
            for (j=0; j<VP_qtdeGenes; j++)
            {
               VP_mapaDist[i][j] = infinito; //Inicia Todos os genes com valor infinito na distância
                                             //ou seja, não tem caminho entre eles
            }
            VP_mapaDist[i][i] = 0.0; //a distância de um gene para ele mesmo é 0
         }
      }
 
      ~TMapaGenes () 
      {
         int i;

         for (i=VP_qtdeGenes-1; i>=0; i--)
         {
            free (VP_mapaDist [i]);
         }

         if (VP_qtdeGenes>0) free (VP_mapaDist);
      }

      void set_distancia(int geneOri, int geneDest, double distancia)
      { 
         //a distância do gene para ele mesmo não pode ser alterada
         //nenum gene pode está fora do indice d tabela
         if ((geneOri!=geneDest)&&(geneOri>=0)&&(geneOri<VP_qtdeGenes)&&(geneDest>=0)&&(geneDest<VP_qtdeGenes))
            VP_mapaDist[geneOri][geneDest] = distancia;         
      }

      double get_distancia(int geneOri, int geneDest)
      {
         //nenum gene pode está fora do indice d tabela
         if ((geneOri<VP_qtdeGenes)&&(geneDest>=0)&&(geneDest<VP_qtdeGenes))
            return VP_mapaDist[geneOri][geneDest];         
         else
            return 0.0;

      }
};


/************************************************************************
   Fim da clase referente ao carregamento da tabela de distâncias
                Início do Algoritmo Genético
*************************************************************************/

/*************************************************************
    Estrutura referente a um gene, representando 1 cidade
Contém 2 ponteriros para o que seria o gene anterior
e posterior.
*************************************************************/
struct TGene
{
   int id;
   struct TGene *prox;
   struct TGene *ant;   
};

/**************************************************
class indivíduo. Uma sequencia de genes encadeados
***************************************************/

class TIndividuo
{
   private:
      int VP_ultOpt;

      TGene **VP_direto;
      TGene *VP_geneIni;
      int VP_qtdeGenes;
      double VP_dist;
      int VP_tipoMutacao;   //0 - 2opt
                            //1 - Simples troca de 2 genes

   /**************************************
    funções privadas (métodos privados)
   ***************************************/
   double calcDistTotTroca (int idG1, int idG2)
   {
      if ((idG1==0)||(idG2==0)||(idG1==idG2)) return VP_dist;

      TGene *tempG;

      double tot=VP_dist;

      //Entre os genes 1 e 2, a soma é do caminho de retorno ou seja, de 2 para 1
      for (tempG=VP_direto[idG1]; tempG->id!=idG2; tempG = tempG->prox)
      {
         tot -= Mapa->get_distancia(tempG->id, tempG->prox->id);
         tot += Mapa->get_distancia(tempG->prox->id, tempG->id); 
      }

      //Arestas que ficaram faltando
      tot -= Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG1]->id);
      tot -= Mapa->get_distancia(VP_direto[idG2]->id, VP_direto[idG2]->prox->id);

      tot += Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG2]->id);
      tot += Mapa->get_distancia(VP_direto[idG1]->id, VP_direto[idG2]->prox->id);

      return tot;
   }

   void mutacaoSimples()
   {
      int g1=0;
      int g2=0;

      g1 = (rand()%(VP_qtdeGenes-1))+1;
      g2 = (rand()%(VP_qtdeGenes-1))+1;
      
      troca(g1, g2);
   }

   void mutacao2OptArtigo()
   {
      int g1=0;
      int g2=0;
      double tot;
      
      TGene *temp;
      TGene *tempG1;
      TGene *tempG2;

      g1 = (rand()%(VP_qtdeGenes-1))+1;
      g2 = (rand()%(VP_qtdeGenes-1))+1;

      tempG1 = VP_direto[g1];
      tempG2 = VP_direto[g2];
      //Fazemdo com que g1 esteja sempre antes de g2
      for (temp=VP_geneIni->prox; temp->id!=0; temp=temp->prox)
      {
         if(g1==g2) break;

         if (tempG1->prox->id==g2) break;
         else if (tempG1->prox->id==0)
         {
            g2 = g1;
            g1 = temp->id;
            break;
         }
         else if (tempG2->prox->id==0) break;
         else if (tempG2->prox->id==g1)
         {
            g2 = g1;
            g1 = temp->id;
            break;
         }
         else if(temp->id==g1) break;
         else if (temp->id==g2)
         {
            g2 = g1;
            g1 = temp->id;
            break;
         }

         tempG1=tempG1->prox;
         tempG2=tempG2->prox;
      }

      tot = calcDistTotTroca (g1, g2);
      if(tot<VP_dist) troca(g1, g2);
   }

   void mutacaoRapida()
   {
      int g1=0;
      int g2=0;
      TGene *tempG;

      g1 = (rand()%(VP_qtdeGenes-1))+1;
      g2 = (rand()%(VP_qtdeGenes-1))+1;

      if(VP_direto[g1]->prox->id==VP_direto[g2]->id)
      {
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ant->id, VP_direto[g1]->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->id, VP_direto[g1]->prox->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->id, VP_direto[g2]->prox->id);

         VP_dist += Mapa->get_distancia(VP_direto[g1]->ant->id, VP_direto[g2]->id);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->id, VP_direto[g1]->id);
         VP_dist += Mapa->get_distancia(VP_direto[g1]->id, VP_direto[g2]->prox->id);

         VP_direto[g1]->ant->prox = VP_direto[g2];
         VP_direto[g2]->ant = VP_direto[g1]->ant;
         VP_direto[g1]->prox = VP_direto[g2]->prox;
         VP_direto[g1]->prox->ant = VP_direto[g1];
         VP_direto[g1]->ant = VP_direto[g2];
         VP_direto[g2]->prox = VP_direto[g1];
      }
      else if (VP_direto[g2]->prox->id==VP_direto[g1]->id)
      {
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ant->id, VP_direto[g2]->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->id, VP_direto[g2]->prox->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->id, VP_direto[g1]->prox->id);

         VP_dist += Mapa->get_distancia(VP_direto[g2]->ant->id, VP_direto[g1]->id);
         VP_dist += Mapa->get_distancia(VP_direto[g1]->id, VP_direto[g2]->id);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->id, VP_direto[g1]->prox->id);

         VP_direto[g2]->ant->prox = VP_direto[g1];
         VP_direto[g1]->ant = VP_direto[g2]->ant;
         VP_direto[g2]->prox = VP_direto[g1]->prox;
         VP_direto[g2]->prox->ant = VP_direto[g2];
         VP_direto[g2]->ant = VP_direto[g1];
         VP_direto[g1]->prox = VP_direto[g2];

      }
      else if (g1!=g2)
      {
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ant->id, VP_direto[g1]->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->id, VP_direto[g1]->prox->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ant->id, VP_direto[g2]->id);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->id, VP_direto[g2]->prox->id);

         VP_dist += Mapa->get_distancia(VP_direto[g1]->ant->id, VP_direto[g2]->id);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->id, VP_direto[g1]->prox->id);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->ant->id, VP_direto[g1]->id);
         VP_dist += Mapa->get_distancia(VP_direto[g1]->id, VP_direto[g2]->prox->id);

         VP_direto[g1]->ant->prox = VP_direto[g2];
         VP_direto[g1]->prox->ant = VP_direto[g2];

         VP_direto[g2]->ant->prox = VP_direto[g1];
         VP_direto[g2]->prox->ant = VP_direto[g1];

         tempG = VP_direto[g1]->prox;
         VP_direto[g1]->prox = VP_direto[g2]->prox;
         VP_direto[g2]->prox = tempG;

         tempG = VP_direto[g1]->ant;
         VP_direto[g1]->ant = VP_direto[g2]->ant;
         VP_direto[g2]->ant = tempG;
      }
   }

   void troca2Opt()
   {
      TGene *g1;
      TGene *g2;
      double temp = VP_dist;

      double tot;

      for (g1= VP_direto[VP_ultOpt]; (g1->prox->prox->id!=0)&&(g1->id!=0); g1=g1->prox)
         for (g2=g1->prox; g2->prox->id!=0; g2=g2->prox)
         {  
            tot = calcDistTotTroca (g1->id, g2->id);
            //printf ("new %f, old %f\n", tot, VP_dist);
            if(tot<VP_dist) 
            {
               troca(g1->id, g2->id);
               g1 =  VP_direto[g2->id];
               //O próximo 2opt não iniciará mais do primeiro 0->prox
               //pois já se sabe que até g2, o caminho já foi melhorado
               //      VP_ultOpt = g2->id; 
               //printf ("\n\nMerece muita atenção - Troquei new %f, old %f\n\n", tot, VP_dist);
               //      return;
            }
         }

    /*  if (temp<VP_dist)
         VP_ultOpt = VP_geneIni->prox->id;
      else
         VP_ultOpt = 0; //Não executará mais o 2opt   */ 
   }
   
   void mutacao3Opt()
   {
      if (VP_ultOpt==0) return;
      if (VP_qtdeGenes<6) return;

      processa3opt();
   }

   //*************
   // mutação 3opt
   //*************
   void processa3opt()
   {
      TGene *g1;
      TGene *g2;
      TGene *g3;

      int melhorTroca = 0;
      for (g1= VP_direto[VP_ultOpt]; g1->prox->prox->prox->prox->prox->id!=0; g1=g1->prox)
         for (g2=g1->prox->prox; g2->prox->prox->prox->id!=0; g2=g2->prox)
         {  
            for (g3=g2->prox->prox; g3->prox->id!=0; g3=g3->prox)
            {
	       // São 4 as possibilidades de troca no 3-Opt. Devemos utilizar
               // a melhor
               melhorTroca = calcMelhorTotTroca3opt (g1->id, g2->id, g3->id );

               if(melhorTroca!=0) 
               {
                  VP_ultOpt = g1->ant->id; 
                  troca3opt(g1->id, g2->id, g3->id, melhorTroca);
                  VP_ultOpt = VP_direto[VP_ultOpt]->prox->id;
                  //O próximo 3opt não iniciará mais do primeiro 0->prox
                  //pois já se sabe que até g1, o caminho já foi melhorado
                  //printf ("Troquei new %f, old %f\n", tot, VP_dist);
                  return;
               }
            }
         }

      //Se saiu do for, é porque não tem o que o 3opt melhorar mais
      VP_ultOpt = 0;
   }
   
   int calcMelhorTotTroca3opt (int idG1, int idG2, int idG3)
   {
      double melhorDist = VP_dist;
      double dist;
      int melhorTroca = 0;

      dist = calcTroca3opt1(idG1, idG2, idG3);
      if (dist<melhorDist) { melhorTroca = 1; melhorDist = dist; }

      dist = calcTroca3opt2(idG1, idG2, idG3);
      if (dist<melhorDist) { melhorTroca = 2; melhorDist = dist; }

      dist = calcTroca3opt3(idG1, idG2, idG3);
      if (dist<melhorDist) { melhorTroca = 3; melhorDist = dist; }

      dist = calcTroca3opt4(idG1, idG2, idG3);
      if (dist<melhorDist) { melhorTroca = 4; melhorDist = dist; }

       return melhorTroca;
   }

   double calcTroca3opt1(int idG1, int idG2, int idG3)
   {
      //Na troca 1, nenhum dos sentidos é alterado, basta subtrair os pesos das 
      //arestas excluidas e adicionar os pesos das arestas incluidas
      double dist;
      
      dist = VP_dist - (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG1]->id)+
                        Mapa->get_distancia(VP_direto[idG2]->ant->id, VP_direto[idG2]->id)+
                        Mapa->get_distancia(VP_direto[idG3]->ant->id, VP_direto[idG3]->id)) +
                       (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG2]->id)+
                        Mapa->get_distancia(VP_direto[idG3]->ant->id, VP_direto[idG1]->id)+
                        Mapa->get_distancia(VP_direto[idG2]->ant->id, VP_direto[idG3]->id));

      return dist;
   }

   double calcTroca3opt2(int idG1, int idG2, int idG3)
   {
      //Na troca 2, os sentidos de g3 ao fim e de g1 ao anterior de g2 não se alteram,
      //porém o sentido de g2 ao anterior de g3 se inverte, passando a ser do anterior 
      // de g3 à g2
      double dist;

      dist = VP_dist - (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG1]->id)+
                        Mapa->get_distancia(VP_direto[idG2]->ant->id, VP_direto[idG2]->id)+
                        Mapa->get_distancia(VP_direto[idG3]->ant->id, VP_direto[idG3]->id)) +
                       (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG3]->ant->id)+
                        Mapa->get_distancia(VP_direto[idG2]->id, VP_direto[idG1]->id)+
                        Mapa->get_distancia(VP_direto[idG2]->ant->id, VP_direto[idG3]->id))+
                        valDistInvert (idG2, VP_direto[idG3]->ant->id);

      return dist;
   }

   double calcTroca3opt3(int idG1, int idG2, int idG3)
   {
      //Na troca 3, os sentidos de g3 ao fim e de g2 ao anterior de g3 não se alteram,
      //porém o sentido de g1 ao anterior de g2 se inverte, passando a ser do anterior 
      // de g2 à g1
      double dist;

      dist = VP_dist - (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG1]->id)+
                        Mapa->get_distancia(VP_direto[idG2]->ant->id, VP_direto[idG2]->id)+
                        Mapa->get_distancia(VP_direto[idG3]->ant->id, VP_direto[idG3]->id)) +
                       (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG2]->id)+
                        Mapa->get_distancia(VP_direto[idG3]->ant->id, VP_direto[idG2]->ant->id)+
                        Mapa->get_distancia(VP_direto[idG1]->id, VP_direto[idG3]->id))+
                        valDistInvert (idG1, VP_direto[idG2]->ant->id);

      return dist;
   }

   double calcTroca3opt4(int idG1, int idG2, int idG3)
   {
      //Na troca 4, aenas o sentido de g3 ao fim não se inverte
      double dist;

      dist = VP_dist - (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG1]->id)+
                        Mapa->get_distancia(VP_direto[idG2]->ant->id, VP_direto[idG2]->id)+
                        Mapa->get_distancia(VP_direto[idG3]->ant->id, VP_direto[idG3]->id)) +
                       (Mapa->get_distancia(VP_direto[idG1]->ant->id, VP_direto[idG2]->ant->id)+
                        Mapa->get_distancia(VP_direto[idG1]->id, VP_direto[idG3]->ant->id)+
                        Mapa->get_distancia(VP_direto[idG2]->id, VP_direto[idG3]->id))+
                        valDistInvert (idG1, VP_direto[idG2]->ant->id) +
                        valDistInvert (idG2, VP_direto[idG3]->ant->id);


      return dist;
   }

   double valDistInvert (int idG1, int idG2)
   {
      double difDist = 0;
      TGene *g;
      for (g = VP_direto[idG1]; g->id!=idG2; g = g->prox )
      {
         difDist += (Mapa->get_distancia(g->prox->id, g->id) - Mapa->get_distancia(g->id, g->prox->id));
      }
      return difDist;
   }

   void troca3opt(int idG1, int idG2, int idG3, int melhorTroca)
   {
      switch (melhorTroca)
      {
         case 1:
            troca3opt1(idG1, idG2, idG3);
            break;

         case 2:
           troca3opt2(idG1, idG2, idG3);
           break;

         case 3:
              troca3opt3(idG1, idG2, idG3);
              break;

         case 4:
              troca3opt4(idG1, idG2, idG3);  
      }
      recalcDist();
   }

   void troca3opt1(int idG1, int idG2, int idG3)
   {
      int aG1 = VP_direto[idG1]->ant->id;
      int aG2 = VP_direto[idG2]->ant->id;
      int aG3 = VP_direto[idG3]->ant->id;

      VP_direto[aG1]->prox = VP_direto[idG2];
      VP_direto[idG2]->ant =  VP_direto[aG1];

      VP_direto[aG2]->prox = VP_direto[idG3];
      VP_direto[idG3]->ant =  VP_direto[aG2];

      VP_direto[aG3]->prox = VP_direto[idG1];
      VP_direto[idG1]->ant =  VP_direto[aG3];
   }

   void troca3opt2(int idG1, int idG2, int idG3)
   {
      int aG1 = VP_direto[idG1]->ant->id;
      int aG2 = VP_direto[idG2]->ant->id;
      int aG3 = VP_direto[idG3]->ant->id;

      inverteSentido (idG2, aG3);

      VP_direto[aG1]->prox = VP_direto[aG3];
      VP_direto[aG3]->ant =  VP_direto[aG1];

      VP_direto[aG2]->prox = VP_direto[idG3];
      VP_direto[idG3]->ant =  VP_direto[aG2];

      VP_direto[idG2]->prox = VP_direto[idG1];
      VP_direto[idG1]->ant =  VP_direto[idG2];
   }

   void troca3opt3(int idG1, int idG2, int idG3)
   {
      int aG1 = VP_direto[idG1]->ant->id;
      int aG2 = VP_direto[idG2]->ant->id;
      int aG3 = VP_direto[idG3]->ant->id;

      inverteSentido (idG1, aG2);

      VP_direto[aG1]->prox = VP_direto[idG2];
      VP_direto[idG2]->ant =  VP_direto[aG1];

      VP_direto[aG3]->prox = VP_direto[aG2];
      VP_direto[aG2]->ant =  VP_direto[aG3];

      VP_direto[idG1]->prox = VP_direto[idG3];
      VP_direto[idG3]->ant =  VP_direto[idG1];
   }

   void troca3opt4(int idG1, int idG2, int idG3)
   {
      int aG1 = VP_direto[idG1]->ant->id;
      int aG2 = VP_direto[idG2]->ant->id;
      int aG3 = VP_direto[idG3]->ant->id;

      inverteSentido (idG1, aG2);
      inverteSentido (idG2, aG3);

      VP_direto[aG1]->prox = VP_direto[aG2];
      VP_direto[aG2]->ant =  VP_direto[aG1];

      VP_direto[idG1]->prox = VP_direto[aG3];
      VP_direto[aG3]->ant =  VP_direto[idG1];

      VP_direto[idG2]->prox = VP_direto[idG3];
      VP_direto[idG3]->ant =  VP_direto[idG2];
   }

   void inverteSentido (int idG1, int idG2)
   {
      TGene *g;  
      TGene *tmpOldProx;
      for (g = VP_direto[idG1]; g->id!=idG2; g = tmpOldProx )
      {
         tmpOldProx = g->prox;
         g->prox = g->ant;
         g->ant = tmpOldProx;
      }
       tmpOldProx = g->prox;
       g->prox = g->ant;
       g->ant = tmpOldProx;
   }

   public:
      TMapaGenes *Mapa;

      TGene *get_ini () { return VP_geneIni; }

      TGene *prox (TGene *gene) { return gene->prox; }
      
      TGene *ant (TGene *gene) { return gene->ant; }

      void set_tipoMutacao (int tipo) { if ((tipo>=0)&&(tipo<=5)) VP_tipoMutacao = tipo; }

      //Cria um novo indivíduo
      void novo ()
      {
         int i;
         VP_tipoMutacao = 0;

         VP_qtdeGenes = Mapa->get_qtdeGenes();
         
         VP_direto = (TGene **) malloc(VP_qtdeGenes*sizeof(TGene *));

         VP_dist = 0;

         VP_geneIni = (TGene *) malloc(sizeof(TGene));
         VP_geneIni->id = 0;
         
         VP_direto[0] = VP_geneIni;
         for (i=1; i<VP_qtdeGenes; i++)
         {
            VP_direto[i] = (TGene *) malloc(sizeof(TGene));
            VP_direto[i]->id = i;
            VP_direto[i-1]->prox = VP_direto[i];
            VP_direto[i]->ant = VP_direto[i-1];
           
            VP_dist += Mapa->get_distancia(i-1, i);
         }

         //Fechando o ciclo
         VP_direto[VP_qtdeGenes-1]->prox = VP_direto[0];
         VP_direto[0]->ant = VP_direto[VP_qtdeGenes-1];
           
         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         VP_dist += Mapa->get_distancia(VP_qtdeGenes-1, 0);
      }

      //distância total
      double get_distancia () { return VP_dist; }

      char *toString ()
      {

         TGene *tempGene;
         char *result;
         char temp[20];

         sprintf(temp, "%d", VP_geneIni->id);
         result = (char *)malloc(7*VP_qtdeGenes*sizeof(char));
         strcpy (result, temp);

         for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
         {
            sprintf (temp, " - %d", tempGene->id);
            strcat (result, temp);            
         }
         return result;
      };     

      //Embaralha os genes de um determinado individuo
      void embaralha ()
      { 

         int i;
         int rd1;
         int rd2;

         int G1A;
         int G1P;
         
         int G2A;
         int G2P;                  

         for (i=0; i<(VP_qtdeGenes/2); i++)
         {
            rd1 = (rand()%(VP_qtdeGenes-1))+1;
            G1A = VP_direto[rd1]->ant->id;
            G1P = VP_direto[rd1]->prox->id;

            rd2 = (rand()%(VP_qtdeGenes-1))+1;
            G2A = VP_direto[rd2]->ant->id;
            G2P = VP_direto[rd2]->prox->id; 

            //Não posso fazer o anterior de rd1 = anterior de rd2, pois ficaria anterior de rd1 = rd1
            if (VP_direto[rd1]->prox->id == VP_direto[rd2]->id)
            {
               VP_direto[rd1]->prox = VP_direto[G2P]; VP_direto[G2P]->ant = VP_direto[rd1];
               VP_direto[rd1]->ant = VP_direto[rd2]; VP_direto[rd2]->prox = VP_direto[rd1];
               VP_direto[rd2]->ant = VP_direto[G1A]; VP_direto[G1A]->prox = VP_direto[rd2];
            }
            //Não posso fazer o posterior de rd1 = posterior de rd2, pois ficaria posterior de rd1 = rd1
            else if (VP_direto[rd1]->ant->id == VP_direto[rd2]->id)
            {
               VP_direto[rd1]->ant = VP_direto[G2A]; VP_direto[G2A]->prox = VP_direto[rd1];
               VP_direto[rd1]->prox = VP_direto[rd2]; VP_direto[rd2]->ant = VP_direto[rd1];
               VP_direto[rd2]->prox = VP_direto[G1P]; VP_direto[G1P]->ant = VP_direto[rd2];
            }
            else
            {
               VP_direto[rd1]->prox = VP_direto[G2P]; VP_direto[G2P]->ant = VP_direto[rd1];
               VP_direto[rd1]->ant = VP_direto[G2A]; VP_direto[G2A]->prox = VP_direto[rd1];
               VP_direto[rd2]->prox = VP_direto[G1P]; VP_direto[G1P]->ant = VP_direto[rd2];
               VP_direto[rd2]->ant = VP_direto[G1A]; VP_direto[G1A]->prox = VP_direto[rd2];
            }

         }
        
         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         recalcDist ();
      } 

      void debug ()
      {
         int i;
         for (i=0; i<VP_qtdeGenes; i++)
         {
           printf ("Gene: %d, prox: %d, ant: %d\n", VP_direto[i]->id, VP_direto[i]->prox->id, VP_direto[i]->ant->id);
         }
      }

      void recalcDist ()
      {
         TGene *tempGene;
         double j;

         VP_dist = 0;
         for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
         {
            VP_dist += Mapa->get_distancia(ant(tempGene)->id, tempGene->id);
         }
         VP_dist += Mapa->get_distancia(ant(VP_geneIni)->id, VP_geneIni->id);
      }

      //Realiza o 2opt
      void mutacao ()
      {
         int rdSort;
         switch (VP_tipoMutacao)
         {
            case 1:
               mutacaoSimples();
               break;

           case 2:
              mutacao2OptArtigo();
              break;

           case 3:
              mutacaoRapida();
              break;
           case 4:
              rdSort = (rand()%(100));
              if      (rdSort < 70) troca2Opt();
              else if (rdSort < 80) mutacaoSimples();
              else if (rdSort < 90) mutacao2OptArtigo();
              else if (rdSort < 100) mutacaoRapida();
              break;
           case 5:
              mutacao3Opt();
              break;
           default:
              troca2Opt();
         }
      }

      //O 2opt
      void troca (int g1, int g2)
      {
         if ((g1==0)||(g2==0)||(g1==g2)) return; //A primeira posição nunca será trocada.

         TGene *tempGene; 
         int vizinho;
      
         bool troca = false;
         int guarda;

         VP_dist = 0;

         if (VP_direto[g1]->prox->id==g2)
         {
            int AG1 = VP_direto[g1]->ant->id;
            int PG2 = VP_direto[g2]->prox->id;

             VP_direto[AG1]->prox = VP_direto[g2]; VP_direto[g2]->ant = VP_direto[AG1];
             VP_direto[g2]->prox  = VP_direto[g1]; VP_direto[g1]->ant = VP_direto[g2];
             VP_direto[g1]->prox  = VP_direto[PG2]; VP_direto[PG2]->ant = VP_direto[g1];

             //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
             VP_ultOpt = VP_geneIni->prox->id;

             recalcDist ();
             return;
         }
  
         else if (VP_direto[g2]->prox->id==g1)
         {
            int AG2 = VP_direto[g2]->ant->id;
            int PG1 = VP_direto[g1]->prox->id;

             VP_direto[AG2]->prox = VP_direto[g1]; VP_direto[g1]->ant = VP_direto[AG2];
             VP_direto[g1]->prox  = VP_direto[g2]; VP_direto[g2]->ant = VP_direto[g1];
             VP_direto[g2]->prox  = VP_direto[PG1]; VP_direto[PG1]->ant = VP_direto[g2];

             //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
             VP_ultOpt = VP_geneIni->prox->id;

             recalcDist ();
             return;
         }

         for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
         {
            //Encontrei um gene para troca
            if ((tempGene->id==g1)||(tempGene->id==g2))
            {
               if (!troca)
               {
                  if(tempGene->id == g1)
                  {
                     guarda = VP_direto[g2]->prox->id;
                     tempGene->ant->prox = VP_direto[g2];
                     VP_direto[g2]->prox = tempGene->ant;
                     tempGene->ant=tempGene->prox;
                     tempGene = VP_direto[g2];
                  }
                  else if(tempGene->id == g2)
                  {
                     guarda = VP_direto[g1]->prox->id;
                     tempGene->ant->prox = VP_direto[g1];
                     VP_direto[g1]->prox = tempGene->ant;
                     tempGene->ant=tempGene->prox;
                     tempGene = VP_direto[g1];
                  }
               }
               else
               {
                  tempGene->prox = VP_direto[guarda];
                  VP_direto[guarda]->ant = tempGene;
               }
                
               troca = !troca;
            }
            
            if(troca)
            {
                vizinho = tempGene->prox->id;
                tempGene->prox = tempGene->ant;
                tempGene->ant = VP_direto[vizinho];
            }

            VP_dist += Mapa->get_distancia(tempGene->ant->id, tempGene->id);
         }

         //Se o indivíduo foi modificado, o 2opt ou 3 opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         VP_dist += Mapa->get_distancia(VP_geneIni->ant->id, VP_geneIni->id);
      }

      TGene *get_gene (int id) { return VP_direto[id]; }

      //Cross Over, transforma o indivíduo em um descendente de outros 2
      void descendente (TIndividuo parceiro1, TIndividuo parceiro2)
      {
         bool *controle;
         TGene *gPar1;
         TGene *gPar2;

         TGene *g1;
         TGene *g2;

         int pivo = (rand()%(VP_qtdeGenes-1))+1;

         //Servirá de controle para os genes já utilzados
         controle = (bool *)calloc(VP_qtdeGenes, sizeof(bool));
         
         gPar1 = parceiro1.get_gene(pivo);
         gPar2 = parceiro2.get_gene(pivo);

         controle[pivo]=true;

         g1 = VP_direto[pivo];
         g2 = VP_direto[pivo];

         bool esq=true;
         bool dir=true;

         int i=1;
          
         //Enquanto não posicionar todos os genes
         while (i<VP_qtdeGenes)
         {
            if (esq)
            {
               if (gPar1->id==0) esq = false;
               else
               {
                  gPar1 = gPar1->ant;
                  if (!controle[gPar1->id])
                  {
                     g1->ant = VP_direto[gPar1->id];
                     VP_direto[gPar1->id]->prox = g1;
                     g1 = VP_direto[gPar1->id];
                     controle[gPar1->id] = true;
                     i++;
                  }
               }
            }
            
            if (dir)
            {
               if (gPar2->prox->id==0) dir = false;
               else
               {
                  gPar2 = gPar2->prox;
                  if (!controle[gPar2->id])
                  {
                     g2->prox = VP_direto[gPar2->id];
                     VP_direto[gPar2->id]->ant = g2;
                     g2 = VP_direto[gPar2->id];
                     controle[gPar2->id] = true;
                     i++;
                  }
               }
            }
            
            if ((!esq)&&(esq==dir))
            {
               gPar1 = gPar1->ant;
               if (!controle[gPar1->id])
               {
                  g1->ant = VP_direto[gPar1->id];
                  VP_direto[gPar1->id]->prox = g1;
                  g1 = VP_direto[gPar1->id];
                  controle[gPar1->id] = true;
                  i++;
               }

               gPar2 = gPar2->prox;
               if (!controle[gPar2->id])
               {
                  g2->prox = VP_direto[gPar2->id];
                  VP_direto[gPar2->id]->ant = g2;
                  g2 = VP_direto[gPar2->id];
                  controle[gPar2->id] = true;
                  i++;
               }
            }
         }

         g1->ant = g2;
         g2->prox = g1;
         free(controle);

         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         recalcDist ();
      }
};

/*********************************************************
classe de População, conterá todos os indivíduos e gerará
as novas gerações
**********************************************************/
class TPopulacao
{
   private:
      double VP_somaDistancias;
      int VP_tamanho;
      TIndividuo **VP_individuos;
      TMapaGenes *VP_mapa;
      int VP_tipoNovaGeracao;       //Como será a substituição da geração antiga pela nova
                                    //0 - os piores
                                    //1 - a menor diferença de distância

       int VP_tipoMutacao;          // 0 - Escolha sempre aleatória
                                    // 1 - O melhor indivíduo nunca sofre mutação
                                    // 2 - Se o melhor sempre sofrerá mutacao

      int VP_tipoMutacaoIndividuo;  //A forma como o individuo sofrerá a mutacao

      //Funções privadas

      /****************************************************
       A função quickSort não está mais sendo usada,
       pois foi substituida pela função qsort da própria
       biblioteca c
       ****************************************************/
      /*
      void quickSort(int ini, int fim)
      {
         TIndividuo *temp;

         //Se tiver apenas 1 elemento, ele está na posição correta
         if (ini>=fim) return;

         //Se tem 2 elentos, ordena os 2, evitando uma rotina mais complexa.
         if (fim==ini+1)
         {
            if (VP_individuos[ini]->get_distancia()>VP_individuos[fim]->get_distancia())
            {
               temp = VP_individuos[ini];
               VP_individuos[ini] = VP_individuos[fim];
               VP_individuos[fim] = temp;
            }          
            return;
         }

         int i = ini;
         int f = fim;

         //chutando um valor para o Pivo
         int valMeio = (fim+ini)/2; 

         double pivo = VP_individuos[valMeio]->get_distancia();

         //Loop do quickSort
         while (f>i)
         {
            while (VP_individuos[i]->get_distancia()<pivo) i++;
            while (VP_individuos[f]->get_distancia()>pivo) f--;

            if (f>=i)
            {
               temp = VP_individuos[i];
               VP_individuos[i] = VP_individuos[f];
               VP_individuos[f] = temp;
               f--;
               i++;
            }
         }

         //Recursividade do quickSort
         quickSort(ini, f);
         quickSort(i, fim);
      }
*/
      //Substitui os piores individuos
      void novaGeracao0(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;

         for (int i=VP_tamanho-qtdePercent; i<VP_tamanho; i++)
         {
            int parceiro1 = (rand()%(VP_tamanho-qtdePercent));
            int parceiro2 = (rand()%(VP_tamanho-qtdePercent));

            // Não permite que o indivíduo cruze com ele mesmo
            if (parceiro1==parceiro2) parceiro2--;
            if (parceiro2<0)          parceiro2 = 1;
 
            VP_somaDistancias -= VP_individuos[i]->get_distancia();
            VP_individuos[i]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);
            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }
      }

      //Substitui os indivíduos com distância mais próxmas
      void novaGeracao1(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;
         int *posindividuo;
         double *difindividuo;

         //Selecionado os individuos co distâncias mais próximas
         //-------------------------------------------------------
         posindividuo = (int *)calloc(qtdePercent, sizeof(int));
         difindividuo = (double *)malloc(qtdePercent*sizeof(double));

         double difDist;
         int indTemp;
         int ind;
         double difTemp;

         for (int i= 1; i<VP_tamanho; i++)
         {
            indTemp = i;
            difDist = VP_individuos[i]->get_distancia() - VP_individuos[i-1]->get_distancia();
            //Se o vetor ainda não foi completamente preenchido, ou se a distancia for menor
            //que a maior que tem na lista, significa que esta diferença e este individuo
            //entrará na lista dos que tem menor diferença de distância
            if ((posindividuo[qtdePercent-1]==0)||(difDist<difindividuo[qtdePercent-1]))
            {
               for (int k=0; k<qtdePercent; k++)
               {
                  if(posindividuo[k]==0)
                  {
                     posindividuo[k] = indTemp;
                     difindividuo[k] = difDist;
                     break;
                  }
                  else if (difindividuo[k]>difDist)
                  {
                     ind = posindividuo[k];
                     difTemp = difindividuo[k];
                     posindividuo[k] = indTemp;
                     difindividuo[k] = difDist;
                     indTemp = ind;
                     difDist = difTemp;
                  }
               }
            }
         }
         // As menores diferenças foram identificadas

         for (int i=0; i<qtdePercent; i++)
         {
            int parceiro1 = (rand()%(VP_tamanho-qtdePercent));
            int parceiro2 = (rand()%(VP_tamanho-qtdePercent));

            // Não permite que o indivíduo cruze com ele mesmo
            // Não permite que um dos parceiros seja o mesmo indivíduo a ser substituido
            while (parceiro1==posindividuo[i])
               parceiro1= parceiro1>=VP_tamanho?0:parceiro1+1;
            while ((parceiro2==posindividuo[i])||(parceiro1==parceiro2)) 
               parceiro2= parceiro2>=VP_tamanho?0:parceiro2+1;


            VP_somaDistancias -= VP_individuos[posindividuo[i]]->get_distancia();
            VP_individuos[posindividuo[i]]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);
            VP_somaDistancias += VP_individuos[posindividuo[i]]->get_distancia();
         }

         free (posindividuo);
         free (difindividuo);
      }


      //Substitui os indivíduos com diferença de distância 0
      void novaGeracao2(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;
         double difDist;
	
         //Substituindo todos que tem diferença de distância 0 para o anterior
         for (int i=1; i<VP_tamanho; i++)
         {
            difDist = VP_individuos[i]->get_distancia() - VP_individuos[i-1]->get_distancia();
            if(difDist == 0)
            {
               VP_somaDistancias -= VP_individuos[i]->get_distancia();
               VP_individuos[i]->set_tipoMutacao(3);
               VP_individuos[i]->mutacao();
               VP_somaDistancias += VP_individuos[i]->get_distancia();
            }
         }

         //Completa o crossover com os piores.
         //Obs.: um dos que tem distância 0 pode tabém ser um dos piores, mas, 
         //      não é preciso se preocupar com isso
         for (int i=VP_tamanho-qtdePercent; i<VP_tamanho; i++)
         {
            int parceiro1 = (rand()%(VP_tamanho-qtdePercent));
            int parceiro2 = (rand()%(VP_tamanho-qtdePercent));

            // Não permite que o indivíduo cruze com ele mesmo
            if (parceiro1==parceiro2) parceiro2--;
            if (parceiro2<0)          parceiro2 = 1;
 
            VP_somaDistancias -= VP_individuos[i]->get_distancia();
            VP_individuos[i]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);
            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }

      }


      //Substitui os indivíduos com diferença de distância 0
      void novaGeracao3(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;
         int qtdeAlt = 0;
         double difDist;
	
         //Substituindo todos que tem diferença de distância 0 para o anterior
         for (int i=1; i<VP_tamanho; i++)
         {
            difDist = VP_individuos[i]->get_distancia() - VP_individuos[i-1]->get_distancia();
            if(difDist == 0)
            {
               int parceiro1=i;
               int parceiro2=i;

               while (parceiro1==i) parceiro1 = (rand()%VP_tamanho); //Não pode fazer parte do crossover
               // Não permite que o indivíduo cruze com ele mesmo
               while ((parceiro2==i)||(parceiro2==parceiro1)) parceiro2 = (rand()%(VP_tamanho-qtdePercent));
 
               VP_somaDistancias -= VP_individuos[i]->get_distancia();
               VP_individuos[i]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);
               VP_somaDistancias += VP_individuos[i]->get_distancia();

	       qtdeAlt++;
            }
         }

         //Completa o crossover com os piores.
         //Obs.: um dos que tem distância 0 pode tabém ser um dos piores, mas, 
         //      não é preciso se preocupar com isso
         for (int i=VP_tamanho-qtdePercent+qtdeAlt; i<VP_tamanho; i++)
         {
            int parceiro1 = (rand()%(VP_tamanho-qtdePercent));
            int parceiro2 = (rand()%(VP_tamanho-qtdePercent));

            // Não permite que o indivíduo cruze com ele mesmo
            if (parceiro1==parceiro2) parceiro2--;
            if (parceiro2<0)          parceiro2 = 1;
 
            VP_somaDistancias -= VP_individuos[i]->get_distancia();
            VP_individuos[i]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);
            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }

      }

      //O melhor individuo sempre sofre mutação
      void mutacao2(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;

         VP_somaDistancias -= VP_individuos[0]->get_distancia();
         VP_individuos[0]->set_tipoMutacao(VP_tipoMutacaoIndividuo);
         VP_individuos[0]->mutacao();
         VP_somaDistancias += VP_individuos[0]->get_distancia();

         for (int i=1; i<qtdePercent; i++)
         {
            int escolha = (rand()%(VP_tamanho-1)+1);

            VP_somaDistancias -= VP_individuos[escolha]->get_distancia();
            VP_individuos[escolha]->set_tipoMutacao(VP_tipoMutacaoIndividuo);
            VP_individuos[escolha]->mutacao();
            VP_somaDistancias += VP_individuos[escolha]->get_distancia();
         }
      }

      //A escolha é completamente aleatória
      void mutacao0(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;

         for (int i=0; i<qtdePercent; i++)
         {
            int escolha = (rand()%(VP_tamanho));

            VP_somaDistancias -= VP_individuos[escolha]->get_distancia();
            VP_individuos[escolha]->set_tipoMutacao(VP_tipoMutacaoIndividuo);
            VP_individuos[escolha]->mutacao();
            VP_somaDistancias += VP_individuos[escolha]->get_distancia();
         }
      }

      //O melhor individuo nunca sofre mutação
      void mutacao1(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;

         for (int i=0; i<qtdePercent; i++)
         {
            int escolha = (rand()%(VP_tamanho-1)+1);

            VP_somaDistancias -= VP_individuos[escolha]->get_distancia();
            VP_individuos[escolha]->set_tipoMutacao(VP_tipoMutacaoIndividuo);
            VP_individuos[escolha]->mutacao();
            VP_somaDistancias += VP_individuos[escolha]->get_distancia();
         }
      }

      //Ignora mutação do indivíduo
      //Melhor indivíduo sofre sempre mutação 2opt
      //Os demais, mutação simples
      void mutacao3(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;

         VP_somaDistancias -= VP_individuos[0]->get_distancia();
         VP_individuos[0]->set_tipoMutacao(0);
         VP_individuos[0]->mutacao();
         VP_somaDistancias += VP_individuos[0]->get_distancia();

         for (int i=1; i<qtdePercent; i++)
         {
            int escolha = (rand()%(VP_tamanho-1)+1);

            VP_somaDistancias -= VP_individuos[escolha]->get_distancia();
            VP_individuos[escolha]->set_tipoMutacao(VP_tipoMutacaoIndividuo);
            VP_individuos[escolha]->mutacao();
            VP_somaDistancias += VP_individuos[escolha]->get_distancia();
         }
      }

      //Ignora mutação do indivíduo
      //Melhor indivíduo sofre sempre mutação 2opt do artigo (aleatória e testa se melhorou)
      //Os demais, mutação simples
      void mutacao4(int percent)
      {
         int qtdePercent = VP_tamanho * percent/100;

         VP_somaDistancias -= VP_individuos[0]->get_distancia();
         VP_individuos[0]->set_tipoMutacao(2);
         VP_individuos[0]->mutacao();
         VP_somaDistancias += VP_individuos[0]->get_distancia();

         for (int i=1; i<qtdePercent; i++)
         {
            int escolha = (rand()%(VP_tamanho-1)+1);

            VP_somaDistancias -= VP_individuos[escolha]->get_distancia();
            VP_individuos[escolha]->set_tipoMutacao(VP_tipoMutacaoIndividuo);
            VP_individuos[escolha]->mutacao();
            VP_somaDistancias += VP_individuos[escolha]->get_distancia();
         }
      }

   public:

      TPopulacao (int tamanho, TMapaGenes *mapa)
      {
         VP_tipoMutacaoIndividuo = 0;
         VP_tipoNovaGeracao = 0;
         VP_tipoMutacao = 0;
         VP_somaDistancias = 0;

         VP_mapa = mapa;

         VP_tamanho = tamanho;

         VP_individuos = (TIndividuo **)malloc(VP_tamanho*sizeof(TIndividuo *));         
         for (int i=0; i<VP_tamanho; i++)
         {
            VP_individuos[i] = new TIndividuo();
            VP_individuos[i]->Mapa = VP_mapa;
            VP_individuos[i]->novo();
            VP_individuos[i]->embaralha();
            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }         
      }

      ~TPopulacao()
      {
         for (int i=VP_tamanho-1; i>=0; i--)
         {
            free (VP_individuos[i]);
         }

         free (VP_individuos);
      }

      void set_tipoNovaGeracao(int tipo)      { if ((tipo>=0)&&(tipo<=3)) VP_tipoNovaGeracao = tipo; }
      void set_tipoMutacao(int tipo)          { if ((tipo>=0)&&(tipo<=4)) VP_tipoMutacao = tipo; }
      void set_tipoMutacaoIndividuo(int tipo) { VP_tipoMutacaoIndividuo = tipo; }

      char *toString ()
      {
         char *result;
         char temp[25];

         result = (char *)malloc(VP_tamanho*(8*VP_mapa->get_qtdeGenes()+25)*sizeof(char));
         strcpy (result, VP_individuos[0]->toString());
         sprintf (temp, " (%f)", VP_individuos[0]->get_distancia());
         strcat (result, temp);           

         for (int i=1; i<VP_tamanho; i++)
         {
            strcat (result, "\n");            
            strcat (result, VP_individuos[i]->toString());   
            sprintf (temp, " (%f)", VP_individuos[i]->get_distancia());
            strcat (result, temp);           
         }
         return result;
      };     

      void ordena () 
      { 
         qsort (VP_individuos, VP_tamanho, sizeof(TIndividuo *), compare);
      //   quickSort (0, VP_tamanho-1); //Meu quickSort
      }  

      void novaGeracao(int percent)
      {
         switch (VP_tipoNovaGeracao)
         {
            case 1:
               novaGeracao1(percent);
               break;

           case 2:
              novaGeracao2(percent);
              break;

           case 3:
              novaGeracao3(percent);
              break;

           default:
              novaGeracao0(percent);
         }

      } 

      void mutacao(int percent)
      {
         switch (VP_tipoMutacao)
         {
            case 1:
               mutacao1(percent);
               break;

           case 2:
              mutacao2(percent);
              break;

           case 3:
              mutacao3(percent);
              break;

           case 4:
              mutacao4(percent);
              break;

           default:
              mutacao0(percent);
         }
      } 

      TIndividuo *get_melhor() { return VP_individuos[0]; }
      TIndividuo *get_pior() { return VP_individuos[VP_tamanho-1]; }
      double distanciaMedia () { return VP_somaDistancias/VP_tamanho; }
};

int compare(const void *x, const void *y)
{
  return ((*(TIndividuo **)x)->get_distancia() - (*(TIndividuo **)y)->get_distancia());
}

/*********************************************************
Classe de confguração, utilizada para acelerar o processo
de testes, estava ruim compilar para cada configuração e
estava ruim digitar todos os parâmetros de confguração
**********************************************************/
class TConfig
{
   private:

   void leInfo(xmlDocPtr doc, xmlNode * a_node)
   {
      xmlNode *cur_node = NULL;
      int gene=0;
      int dist;
      xmlChar *key;
      int val;

      for (cur_node = a_node; cur_node; cur_node = cur_node->next) 
      {
         if (cur_node->type == XML_ELEMENT_NODE) 
         {
            key = xmlNodeListGetString(doc, cur_node->xmlChildrenNode, 1);
            val = atoi((char *)key);
            xmlFree(key); 

            if (!xmlStrcmp(cur_node->name, (xmlChar *)"tamanhoPopulacao")) tamPopulacao = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"numGeracoes")) maxGeracao = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentManipulacao")) percentManipulacao = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentMutacao")) percentMutacao = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"tipoMutacaoPopulacao")) tipoMutacaoPopulacao = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"tipoMutacaoIndividuo")) tipoMutacaoIndividuo = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"tipoNovaGeracao")) tipoNovaGeracao = val;
            else if (!xmlStrcmp(cur_node->name, (xmlChar *)"printParcial")) printParcial = val;
         }

         leInfo(doc, cur_node->children);
      }
   }

   public:
      /*******************************************
           Parâmetros de teste do sistema
      ********************************************/
      int tamPopulacao;
      int maxGeracao;
      int percentManipulacao;   //percentual de manipulação do indivíduo. (exclusão / cruzamento)
      int percentMutacao;       //percentual de mutação
      int tipoMutacaoPopulacao; //forma de mutação a ser executada na população, considerando o melhor elemento
      int tipoMutacaoIndividuo; //tipo de mutação do indivíduo, 2Opt ou simples
      int tipoNovaGeracao;      //forma de geração da nova geração
      int printParcial;         //Se inprime informações intermediárias

      /*******************************************************
           Os valores padrões são os utilizados no artigo
      ********************************************************/
      TConfig()
      {
         tamPopulacao = 200;
         maxGeracao = 300;
         percentManipulacao = 30;
         percentMutacao = 20;
         tipoMutacaoPopulacao = 0;
         tipoMutacaoIndividuo = 0;
      }

      void carregaDoArquivo(char *nomeArquivo)
      {
         int numGenes;
         xmlDoc *doc = NULL;
         xmlNode *root_element = NULL;

         // Lendo o arquivo 
         doc = xmlReadFile(nomeArquivo, NULL, 0);

         if (doc == NULL) 
         {
            printf("Erro ao carregar o arquivo %s\n", nomeArquivo);
            return;
         }

         // Obtendo o elemento root
         root_element = xmlDocGetRootElement(doc);
         
         //preenchendo a tabela com os valores da distáncia
         leInfo(doc, root_element);

         //liberando documento
         xmlFreeDoc(doc);
         // liberando as variaveis lobais
         xmlCleanupParser();

      }
     

};

/*************************************************************************
 Principal, o arquivo xml a ser carregado deve ser passado como parâmetro

 Linha de execução
 1) executa o programa com as confgurações padrões
    tsp 

 2) executa o programa com interação sobre os valores de configuração
    tsp -c  

 3) executa o programa lendo um arquivo de configuração em xml
    tsp -c <arquivo de confguração> 
**************************************************************************/
int main(int argc, char **argv)
{
   TMapaGenes mapa ;
   TIndividuo *melhor;
   TIndividuo *pior; 
   TConfig config; 

   if (argc < 2)
      return(1);

   LIBXML_TEST_VERSION

   mapa.carregaDoArquivo (argv[1]);

   if (argc>2)
   {
      if ((argv[2][0] == '-')&&(argv[2][1] == 'c'))
      {
         if (argc > 3)
         {
            config.carregaDoArquivo(argv[3]);
         }
         else
         {
            printf("Mostra informações parciais: "); scanf("%d", &config.printParcial);
            printf("Tamanho da população: "); scanf("%d", &config.tamPopulacao);
            printf("Número máximo de gerações: "); scanf("%d", &config.maxGeracao);
            printf("Percentual de manipuação da população: "); scanf("%d", &config.percentManipulacao);
            printf("Percentual de mutação: "); scanf("%d", &config.percentMutacao);
            printf("Forma de mutação da população: "); scanf("%d", &config.tipoMutacaoPopulacao);
            printf("Tipo de mutação do indivíduo: "); scanf("%d", &config.tipoMutacaoIndividuo);
            printf("Forma como a nova geração será gerada: "); scanf("%d", &config.tipoNovaGeracao);
         }
      }
      else
         return 1;
   }

    time_t sysTime1, sysTime2;
    time(&sysTime1);
    TPopulacao populacao (config.tamPopulacao, &mapa);
    populacao.set_tipoMutacao(config.tipoMutacaoPopulacao);
    populacao.set_tipoNovaGeracao (config.tipoNovaGeracao);
    populacao.mutacao (config.percentMutacao);
    populacao.set_tipoMutacaoIndividuo(config.tipoMutacaoIndividuo);
    populacao.ordena();

    for(int i=0; i<config.maxGeracao; i++)
    {
       populacao.novaGeracao (config.percentManipulacao);
       populacao.mutacao (config.percentMutacao);
       populacao.ordena();
       if(config.printParcial==1)
       {
          melhor = populacao.get_melhor();
          pior = populacao.get_pior();
          time(&sysTime2);
          printf ("Geração %d. Melhor : %f - Pior: %f - Media: %f (%f)\n",i, melhor->get_distancia(), pior->get_distancia(), populacao.distanciaMedia(), difftime(sysTime2, sysTime1));
       }
    }

    time(&sysTime2);
    printf ("Temmpo de execução %f\n", difftime(sysTime2, sysTime1));
    
    melhor = populacao.get_melhor();
    printf ("%s (%f)\n", melhor->toString(), melhor->get_distancia());
     
/*
    TIndividuo ind;
    TIndividuo ind2;
    TIndividuo ind3;

    ind.Mapa = &mapa;
    ind.novo();
    ind.embaralha();

    ind2.Mapa = &mapa;
    ind2.novo();
    ind2.embaralha();

    printf ("ind1 %s\n", ind.toString());
    printf("%f\n", ind.get_distancia());
    ind.debug();

    printf ("\n");

    printf ("ind2 %s\n", ind2.toString());
    printf("%f\n", ind2.get_distancia());
    ind2.debug();

    printf ("\n");
    ind3.Mapa = &mapa;
    ind3.novo();
    ind3.descendente(ind, ind2);

    printf ("ind3 %s\n", ind3.toString());
    printf("%f\n", ind3.get_distancia());
    ind3.debug();

/*

    ind.mutacao();

    printf ("%s\n", ind.toString());
    ind.debug();
    printf("%f\n", ind.get_distancia());


    printf("%f\n", mapa.get_distancia(0,2));
    printf("%f\n", mapa.get_distancia(2,2));
    printf("%f\n", mapa.get_distancia(1,2));
    printf("%f\n", mapa.get_distancia(2,1));
    printf("%f\n", mapa.get_distancia(50,1));
    printf("%f\n", mapa.get_distancia(50,49));
*/
    return 0;
}
#else
int main(void) {
    fprintf(stderr, "Suporte a árvore xml não compilado\n");
    exit(1);
}
#endif
