## Gromacs file: .gro
## This is for personal use..
## It cleans up residue name, residue name number
## very useful for detecting wrong name order and missin atom in residue


def name_organazier(file_name):

  messy_file = file_name.copy()

  res_name_num = messy_file['res_name'].str.extract('([0-9]+)?')
  res_name = messy_file['res_name'].str.extract('([a-zA-Z]+)')

  messy_file['res_name'] = res_name
  messy_file.insert(0, "res_name_num", res_name_num)
  messy_file['res_name_num'] = pd.to_numeric(messy_file['res_name_num'], errors='coerce')

  messy_name_list = messy_file['res_name'].unique().tolist()

  clean_file = []

  for i in range(0,len(messy_name_list)):

    messy_name_tempo = messy_file[messy_file['res_name'] == messy_name_list[i]]

    clean_file.append(messy_name_tempo)

  n = 1

  clean_file = pd.concat(clean_file, ignore_index=True)

  final_file = []

  for j in range(0,len(messy_name_list)):

    messy_name_tempo = clean_file[clean_file['res_name'] == messy_name_list[j]]

    messy_name_num = list(np.repeat(range(n,n+int(len(messy_name_tempo['res_name_num'])/len(set(messy_name_tempo['res_id'])))),
                                    len(set(messy_name_tempo['res_id']))))

    messy_name_tempo['res_name_num'] = messy_name_num

    n += len(set(messy_name_tempo['res_name_num']))

    final_file.append(messy_name_tempo)

  final_file = pd.concat(final_file, ignore_index=True)

  final_file['res_num'] = range(1,len(final_file)+1)

  # Box
  box_x = max(final_file['x'])
  box_y = max(final_file['y'])
  box_z = max(final_file['z'])

  box_vectors = [box_x, box_y, box_z]

  # GRO file format
  format_string = "{:5}{:<5}{:>5}{:5}{:8.3}{:8.3}{:8.3}"
  format_string_num = "{:5}{:<5}{:>5}{:5}{:8.3f}{:8.3f}{:8.3f}"

  # Open a file and write the header
  with open('output.gro', 'w') as file:
    header = format_string.format('res_num', 'res_name', 'res_id', 'res_num', 'x', 'y', 'z')
    file.write('frame t=1.000 '+ '\n')
    file.write('{:5}'.format(max(final_file['res_num']))+ '\n')

    # Write each row using the format string
    for index, row in final_file.iterrows():
        row_string = format_string_num.format(row['res_name_num'], str(row['res_name']), str(row['res_id']), row['res_num'],
                                          row['x'], row['y'], row['z'])
        file.write(row_string + '\n')

    file.write('  {:8.5f}   {:8.5f}   {:8.5f}\n'.format(*box_vectors))

  return final_file
