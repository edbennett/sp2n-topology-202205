#!/usr/bin/env bash

set -eu

SCALES_FILE=$1
T0_FILE=$2
W0_FILE=$3

cat > ${T0_FILE} <<'EOF'
\begin{tabular}{|c|c|c|c|}
\hline
\hline
$\beta$ & $\mathcal{E}_0$ & $\sqrt{t_0}/a$ (pl.) & $\sqrt{t_0}/a$ (cl.) \\
\hline
EOF

cat > ${W0_FILE} <<'EOF'
\begin{tabular}{|c|c|c|c|}
\hline
\hline
$\beta$ & $\mathcal{W}_0$ & $w_0/a$ (pl.) & $w_0/a$ (cl.) \\
\hline
EOF

awk '{print "$",$3,"$&$",$4,"$&$",$5,"$&$",$6,"$ \\\\"}' ${SCALES_FILE} >> ${T0_FILE}
awk '{print "$",$3,"$&$",$7,"$&$",$8,"$&$",$9,"$ \\\\"}' ${SCALES_FILE} >> ${W0_FILE}

cat >> ${T0_FILE} <<EOF
\hline
\end{tabular}
EOF

cat >> ${W0_FILE} <<EOF
\hline
\end{tabular}
EOF
