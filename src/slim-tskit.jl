function toslim(model::MainlandIslandModel, genetic_map, ngen; 
        tsout=nothing, print_freqs=false, log=0, N1=1)
    @warn "Assuming constant recombination and mutation rates!"
    @unpack arch, m, N = model
    str  = "initialize() {"
    if !isnothing(tsout)
        str *= "\ninitializeTreeSeq();\n"
    end
    imap = indexmap(arch.loci)
    r    = Barriers.recrate(genetic_map, 1, 2)
    u    = arch[1].u
    ys   = [genetic_map[arch.xs[i]] for i=1:length(arch)] 
    map(collect(imap)) do (k,i)
        str *= """
        \tinitializeMutationRate($u);
        \tinitializeMutationType("m$i", $(k.h), "f", $(-k.s));
        \tinitializeGenomicElementType("g$i", m$i, 1.0);
        """
    end
    insert = map(1:length(arch)) do i
        locus = arch[i]
        ix  = imap[locus]
        pos = genetic_map[arch.xs[i]] - 1
        str *= """
        \tinitializeGenomicElement(g$ix, $pos, $pos);
        """
        "p1.genomes.addNewDrawnMutation(m$ix, $pos);"
    end
    str *= """
    initializeRecombinationRate($r, $(physlength(genetic_map)));
}
1 early() {
    sim.addSubpop("p1", $N1);
    sim.addSubpop("p2", $N);
    p1.setMigrationRates(p2, 0.0);  // p2->p1  no migration
    p2.setMigrationRates(p1, $m); // p1->p2
}
1 early() {
    $(join(insert, "\n\t"))
}
// mutationEffect(m2, p1) { return 0.0; }  // This is how one would model divergent sel
mutation(NULL, p1) { return F; }  // see SLiM 4.0.1 manual p. 725
1:$ngen late() { 
"""
    if !isnothing(tsout)
        str *= "sim.treeSeqOutput(\"$tsout\");\n"
    end
    if print_freqs
        str *= "print(sim.mutationFrequencies(p2, sim.mutationsOfType(m1)));\n"
    end
    str *= "\n}\n"
    if log > 0
        str *= "1:$ngen early() { if (sim.cycle % $log == 1){ print(sim.cycle); }}"
    end
    return str 
end

function ts_get_mutations(ts)
    map(ts.variants()) do var
        gt = var.genotypes.tolist()
        var.site.position, classify(gt)
    end
end

