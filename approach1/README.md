# Project Approach 1: Divide and Conquer pplacer

\begin{algorithm}[H]
\SetKwFor{ParallelFor}{parallel for}{do}{endfor}
\SetAlgoLined
\KwResult{$T'$, tree T with query sequence $q$ added}
\KwIn{Tree $T$ on $N$ sequences, the MSA of $N+1$ sequences, and query sequence $q$}
 // $\operatorname{centroidDecomposition}$ decomposes a tree into roughly equal size, disjoint parts.\;
 // $\operatorname{modifyTree}$ adds sequence to a tree based on the sequence's location in the a subtree with the sequence added\;
 $\{T_1,\dots,T_n\} \leftarrow \operatorname{centroidDecomposition}(T)$\;
 $\{S_1, \dots, S_n\} \leftarrow 0$ // Score for each tree\;
 \ParallelFor{$i=1,\dots,n$}{
  // Place query sequence $q$ into the subtree\;
  $T'_i \leftarrow \operatorname{pplacer}(T_i, q)$\;
  // Add the location of the query sequence $q$ to a copy of $T$\;
  $T_{q_{i}} \leftarrow \operatorname{modifyTree}(T, T'_i, q)$\;
  // $\operatorname{RAxMLScorer}$ runs RAxML in fixed tree mode.\;
  // The output score is the maximum likelihood found on the tree.\;
  $S_i \leftarrow \operatorname{RAxMLScorer}( T_{q_{i}})$\;
 }
 // Do a maxLoc reduction for the tree\;
 $bestTreeIndex \leftarrow \operatorname{argmax}_{i} (S_1,\dots,S_n)$\;
 return $T'_{q_{bestTreeIndex}}$\;
 \caption{divide-and-conquer pplacer}
\end{algorithm}

