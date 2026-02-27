for i in $(find . -type f -name "*.nim");
do
   if [[ -f "${i%.nim}" ]]; then
     echo removing "${i%.nim}"
     rm "${i%.nim}"
   fi
done
