#!/usr/bin/env python
"""
PathExplainerTorch - local copy from suinleelab/path_explain.
Source: https://github.com/suinleelab/path_explain/blob/master/path_explain/explainers/path_explainer_torch.py
"""
import functools
import operator
import torch
from torch.autograd import grad
import numpy as np
from tqdm import tqdm


def gather_nd(params, indices):
    """
    Args:
        params: Tensor to index
        indices: k-dimension tensor of integers.
    Returns:
        output: 1-dimensional tensor of elements of ``params``, where
            output[i] = params[i][indices[i]]
    """
    max_value = functools.reduce(operator.mul, list(params.size())) - 1
    indices = indices.t().long()
    ndim = indices.size(0)
    idx = torch.zeros_like(indices[0]).long()
    m = 1

    for i in range(ndim)[::-1]:
        idx += indices[i] * m
        m *= params.size(i)

    idx[idx < 0] = 0
    idx[idx > max_value] = 0
    return torch.take(params, idx)


class PathExplainerTorch(object):
    def __init__(self, model):
        self.model = model
        return

    def _get_ref_tensor(self, baseline, batch_size, num_samples):
        number_to_draw = num_samples * batch_size
        replace = baseline.shape[0] < number_to_draw
        sample_indices = np.random.choice(
            baseline.shape[0], size=number_to_draw, replace=replace
        )
        ref_tensor = baseline[sample_indices, :]
        return ref_tensor

    def _get_samples_input(self, input_tensor, baseline, num_samples, use_expectation):
        input_dims = list(input_tensor.size())[1:]
        num_input_dims = len(input_dims)
        batch_size = input_tensor.size()[0]

        if use_expectation:
            reference_tensor = self._get_ref_tensor(baseline, batch_size, num_samples)
            shape = reference_tensor.shape
            reference_tensor = reference_tensor.view(
                batch_size, num_samples, *(shape[1:])
            )
            t_tensor = torch.FloatTensor(batch_size, num_samples).uniform_(0, 1).to(
                reference_tensor.device
            )
            shape = [batch_size, num_samples] + [1] * num_input_dims
            interp_coef = t_tensor.view(*shape)
            end_point_ref = (1.0 - interp_coef) * reference_tensor
            input_expand_mult = input_tensor.unsqueeze(1)
            end_point_input = interp_coef * input_expand_mult
            samples_input = end_point_input + end_point_ref
        else:
            batch_size = input_tensor.size()[0]
            input_expand = input_tensor.unsqueeze(1)
            reps = np.ones(len(baseline.shape)).astype(int)
            reps[0] = batch_size
            reference_tensor = baseline.repeat(list(reps)).unsqueeze(1)
            scaled_inputs = [
                reference_tensor
                + (float(i) / (num_samples - 1)) * (input_expand - reference_tensor)
                for i in range(0, num_samples)
            ]
            samples_input = torch.cat(scaled_inputs, dim=1)

        samples_delta = self._get_samples_delta(input_tensor, reference_tensor)
        samples_delta = samples_delta.to(samples_input.device)
        return samples_input, samples_delta

    def _get_samples_delta(self, input_tensor, reference_tensor):
        input_expand_mult = input_tensor.unsqueeze(1)
        sd = input_expand_mult - reference_tensor
        return sd

    def _get_grads(self, samples_input, output_indices=None):
        grad_tensor = torch.zeros(samples_input.shape).float().to(samples_input.device)
        k_ = samples_input.shape[1]

        for i in range(k_):
            particular_slice = samples_input[:, i]
            batch_output = self.model(particular_slice)
            if batch_output.size(1) > 1:
                sample_indices = torch.arange(0, batch_output.size(0)).to(
                    samples_input.device
                )
                indices_tensor = torch.cat(
                    [sample_indices.unsqueeze(1), output_indices.unsqueeze(1)], dim=1
                )
                batch_output = gather_nd(batch_output, indices_tensor)

            model_grads = grad(
                outputs=batch_output,
                inputs=particular_slice,
                grad_outputs=torch.ones_like(batch_output).to(samples_input.device),
                create_graph=True,
            )
            grad_tensor[:, i, :] = model_grads[0]
        return grad_tensor

    def attributions(
        self, input_tensor, baseline, num_samples=50, use_expectation=True, output_indices=None
    ):
        equal_dims = baseline.shape[1:] == input_tensor.shape[1:]
        almost_equal_dims = baseline.shape == input_tensor.shape[1:]

        dev = input_tensor.device
        baseline = baseline.to(dev)
        input_tensor.requires_grad_ = True

        if use_expectation and not equal_dims:
            raise ValueError(
                "baseline should be shape (num_refs, ...) "
                "where ... indicates the dimensionality of the input"
            )

        if not use_expectation and baseline.shape[0] != 1:
            if almost_equal_dims:
                baseline = baseline.unsqueeze(0)
            else:
                raise ValueError(
                    "baseline should be shape (...) "
                    "where ... indicates the dimensionality of the input"
                )

        samples_input, samples_delta = self._get_samples_input(
            input_tensor, baseline, num_samples, use_expectation
        )
        grad_tensor = self._get_grads(samples_input, output_indices)
        mult_grads = samples_delta * grad_tensor
        attributions = mult_grads.mean(1)
        return attributions
