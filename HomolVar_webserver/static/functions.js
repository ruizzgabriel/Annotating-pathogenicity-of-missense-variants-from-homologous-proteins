/* user_parameters.html */
function showAminoAcidList(input) {
    const aminoAcidList = document.getElementById('aminoAcidList');
    aminoAcidList.style.display = 'block';
    const rect = input.getBoundingClientRect();
    aminoAcidList.style.top = rect.bottom + window.scrollY + 'px';
    aminoAcidList.style.left = rect.left + window.scrollX + 'px';

    // Hide list when clicking outside
    document.addEventListener('click', function(event) {
        if (!input.contains(event.target) && !aminoAcidList.contains(event.target)) {
            aminoAcidList.style.display = 'none';
        }
    });
}

function selectAminoAcid(aminoAcid) {
    document.getElementById('input1').value = aminoAcid;
    document.getElementById('input3').value = aminoAcid;  // Set both input1 and input3 to the selected amino acid
    document.getElementById('aminoAcidList').style.display = 'none';
}

