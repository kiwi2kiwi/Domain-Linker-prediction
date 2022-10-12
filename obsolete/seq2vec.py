input_file = "C:\\Users\\yanni\\Downloads\\dssp-sequence.fasta"
stream = open(input_file,"r")
fasta_dict = {}
id = ""
for line in stream:
    if line.startswith(">"):
        id = line[1:].strip()
        if id not in set:
            set.add(id)
            fasta_dict[id] = fasta_protein(id)
            fasta_dict[id].range = range




plt.close()
predictions, labels = test_output(training_dataset_list[4], network_es, device)
make_evaluation(predictions, labels, options.remote, ("default" + "default cnn"))
df = pd.DataFrame(train_losses, columns=['loss'])
df['x'] = [i for sl in [[i] * int(len(train_losses) / len(test_losses)) for i in range(len(test_losses))]
           for i
           in sl]
df['type'] = 'train'
df2 = pd.DataFrame((enumerate(test_losses)), columns=['x', 'loss'])
df2['type'] = 'validation_set'
df = pd.concat((df, df2)).reset_index()
fig = sns.relplot(data=df, x='x', y='loss',
                  row=False, row_order=['train', 'validation_set'],
                  kind='line', lw=2.5, aspect=2,
                  hue='type', palette='mako')
fig.set(xlim=(0, max(df.x)), ylim=(0, 1), xlabel='Epoch',
        ylabel="1 - "+plot_stats, xticks=list(range(max(df.x))))
fig.set_titles(row_template='{row_name}')
[l.set_linewidth(2.5) for l in fig.legend.legendHandles]
fig.tight_layout()
sns.despine()
plt.show()
