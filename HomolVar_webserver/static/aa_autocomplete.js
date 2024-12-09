$(document).ready(function() {
    // Arrays for initial and final amino acids
    const initialAminoAcids = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'];
    const finalAminoAcids = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', '*/Ter/fs', '+'];

    // Populate the datalists with options
    const populateDatalist = (datalistId, aminoAcids) => {
        const datalist = document.getElementById(datalistId);
        aminoAcids.forEach(acid => {
            const option = document.createElement('option');
            option.value = acid;
            datalist.appendChild(option);
        });
    };

    // Populate initial amino acids datalists
    populateDatalist('initial-amino-acids-1', initialAminoAcids);
    populateDatalist('initial-amino-acids-2', initialAminoAcids);

    // Populate final amino acids datalists
    populateDatalist('final-amino-acids-1', finalAminoAcids);
    populateDatalist('final-amino-acids-2', finalAminoAcids);

    // Autocomplete behavior for initial amino acids (i_aa1 and i_aa2)
    $('input[name="i_aa1"], input[name="i_aa2"]').on('input', function() {
        const inputValue = $(this).val().toLowerCase();
        const filteredAminoAcids = initialAminoAcids.filter(acid =>
            acid.toLowerCase().startsWith(inputValue)
        );
        const datalistId = $(this).attr('list');
        const datalist = document.getElementById(datalistId);
        datalist.innerHTML = '';
        filteredAminoAcids.forEach(acid => {
            const option = document.createElement('option');
            option.value = acid;
            datalist.appendChild(option);
        });
    });

    // Autocomplete behavior for final amino acids (f_aa1 and f_aa2)
    $('input[name="f_aa1"], input[name="f_aa2"]').on('input', function() {
        const inputValue = $(this).val().toLowerCase();
        const filteredAminoAcids = finalAminoAcids.filter(acid =>
            acid.toLowerCase().startsWith(inputValue)
        );
        const datalistId = $(this).attr('list');
        const datalist = document.getElementById(datalistId);
        datalist.innerHTML = '';
        filteredAminoAcids.forEach(acid => {
            const option = document.createElement('option');
            option.value = acid;
            datalist.appendChild(option);
        });
    });
});
