#!/usr/bin/env bash

set -eu

cat <<'EOF'
\begin{tabular}{|c|cc|cc|}
\hline
\hline
$N_c$ & $\chi_L t_0^2(a=0) & \tilde{\mathcal{X}}^2 & $\chi_L t_0^2(a=0) & \tilde{\mathcal{X}}^2 \\
 & $c_e = 0.225$ & & c_e = 0.5 & \\
\hline
EOF

paste -d' ' <(printf "\$2\$\n\$4\$\n\$6\$\n\$8\$\n") $1 $2 <(printf '\\\\\n\\\\\n\\\\\n\\\\\n')

cat <<'EOF'
\hline
\hline
$N_c$ & $\chi_L w_0^4(a=0) & \tilde{\mathcal{X}}^2 & $\chi_L w_0^4(a=0) & \tilde{\mathcal{X}}^2 \\
 & $c_w = 0.225$ & & c_w = 0.5 & \\
\hline
EOF

paste -d' ' <(printf "\$2\$\n\$4\$\n\$6\$\n\$8\$\n") $3 $4 <(printf '\\\\\n\\\\\n\\\\\n\\\\\n')

cat <<'EOF'
\hline
\end{tabular}
EOF
