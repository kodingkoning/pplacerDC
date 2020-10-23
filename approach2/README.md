
# Approach 2: APPLES with pplacer

\begin{algorithm}[H]
\SetAlgoLined
\KwResult{$T'$, tree T with query sequence $q$ added}
\KwIn{Tree $T$ on $N$ sequences, the MSA of $N+1$ sequences, and query sequence $q$}
 // $\operatorname{centroidDecomposition}$ decomposes a tree into roughly equal size, disjoint parts.\;
 // $\operatorname{modifyTree}$ adds sequence to a tree based on the sequence's location in the a subtree with the sequence added\;
 // $\operatorname{getRegion}$ finds a subtree containing the sequence with a maximum number of sequences\;
  // Run APPLES\;
  $T'_{APPLES} \leftarrow \operatorname{APPLES}(T, q)$\;
  // Identify region of T where q was placed with fewer than 1000 sequences\;
  $T_{qregion} \leftarrow \operatorname{getRegion}(T'_{APPLES}, q, 1000)$\;
  // Run pplacer on the area of T around the placement of q\;
  $T'_{pplacer} \leftarrow \operatorname{pplacer}(T_{qregion})$\;
  return $T'_{pplacer}$\;
\caption{APPLES with pplacer}
\end{algorithm}

