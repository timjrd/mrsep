# Diversité des souches de pathogènes 
## Cedric Chauve (cedric.chauve@sfu.ca)
## 21/03/2019 

Ce repo github contient le travail développé au LaBRI en 2019 pour un
projet visant à estimer la diversité des souches d'un pathogène
présentes dans un jeu de données de lectures d'ADN (1).

Ce repo contiendra tout ce qui permet de suivre nos travaux,
d'analyser et de reproduire nos résultats.

La structure du repo est la suivante.

- Le répertoire src doit contenir les sources de tous les programes
  que nous allons utiliser, incluant les programmes externes et les
  programmes que nous fournirons comme "produit" final. Dans l'idéal,
  ce répertoire ne doit PAS contenir les versions intermediaires des
  programmes que nous allons dévélopper.

- Le répertoire bin contient les exécutables des programmes externes
  dans le cas où nous ne les avons pas compilés nous-même.

- Le répertoire dev est utilisé pour développer des prototypes de nos
  programmes, effectuer de petits tests. Le contenu de de répertoire
  n'a pas vocation à être distribué, ni conservé.

- Le prépertoire doc contient la documentation finale associée à notre
  projet.
  - Le sous-répertoire doc/notebooks contient des notebook Jupyter
    (https://github.com/jupyterlab/jupyterlab). 
  - Le sous-repertoire doc/papers contient les articles (références)
    reliés à notre projet.
  - Le sous-répertoire doc/notes contient les notes et documents
    préliminaires/intermédiaires.

- Le répertoire data contient les données que nous allons utiliser et
  qui ne seront pas modifiées au cours du projet.

- Le répertoire exp contient toutes les expériences que nous allons
  mener. Il est crucial de suivre des règles très strictes sur la
  structure de ce répertoire si nous voulons être capable de suivre
  nos progrès et de comprendre ce que nous faisons. Cette structure
  est expliquée dans le fichier exp/README_structure.

Chaque répertoire a un fichier README qui décrit son contenu et permet
de naviguer dans le projet. Iol est IMPERATIF que, dès qu'un nouveau
répertoire est créé, le fichier README de son parent soit mis à jour
et un fichier README soit créé dans le nouveau répertoire.

(1) A l'exception de ce fichier README.md, l'ensemble de la
documentation est en anglais.